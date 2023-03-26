

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>

#include "helpers/debug.hpp"
#include "helpers/endswith.hpp"
#include "io/jsonhandler.hpp"
#include "model/model.hpp"

// Setting scalar type to double
typedef double T;

int main(int argc, char *argv[]) {
  DEBUG_MSG("Starting svZeroDCalibrator");

  // Get input and output file name
  if (argc != 3) {
    std::cout << "Usage: calibrator path/to/config.json path/to/output.csv"
              << std::endl;
    exit(1);
  }
  std::string input_file = argv[1];
  std::string output_file = argv[2];

  std::ifstream input_file_stream(input_file);
  std::stringstream buffer;
  buffer << input_file_stream.rdbuf();
  std::string config = buffer.str();

  // Load model and configuration
  DEBUG_MSG("Read configuration");
  auto handler = IO::JsonHandler(config);

  auto model = new MODEL::Model<T>();
  std::vector<std::tuple<std::string_view, std::string_view>> connections;

  // Create vessels
  DEBUG_MSG("Load vessels");
  std::map<std::int64_t, std::string_view> vessel_id_map;
  auto vessels = handler["vessels"];
  int param_counter = 0;
  for (size_t i = 0; i < vessels.length(); i++) {
    auto vessel_config = vessels[i];
    std::string_view vessel_name = vessel_config.get_string("vessel_name");
    vessel_id_map.insert({vessel_config.get_int("vessel_id"), vessel_name});
    model->add_block(MODEL::BlockType::BLOODVESSEL,
                     {param_counter, param_counter + 1, param_counter + 2},
                     vessel_name);
    DEBUG_MSG("Created vessel " << vessel_name);
    param_counter += 3;
  }

  // Create junctions
  auto junctions = handler["junctions"];
  for (size_t i = 0; i < junctions.length(); i++) {
    auto junction_config = junctions[i];
    auto j_type = junction_config.get_string("junction_type");
    auto junction_name = junction_config.get_string("junction_name");
    if ((j_type == "NORMAL_JUNCTION") || (j_type == "internal_junction")) {
      model->add_block(MODEL::BlockType::JUNCTION, {}, junction_name);
    } else if (j_type == "BloodVesselJunction") {
      auto outlet_vessels = junction_config.get_int_array("outlet_vessels");
      int num_outlets = outlet_vessels.size();
      std::vector<int> param_ids;
      for (size_t i = 0; i < num_outlets * 3; i++) {
        param_ids.push_back(param_counter++);
      }
      model->add_block(MODEL::BlockType::BLOODVESSELJUNCTION, param_ids,
                       junction_name);
    } else {
      throw std::invalid_argument("Unknown junction type");
    }
    // Check for connections to inlet and outlet vessels and append to
    // connections list
    for (auto vessel_id : junction_config.get_int_array("inlet_vessels")) {
      connections.push_back({vessel_id_map[vessel_id], junction_name});
    }
    for (auto vessel_id : junction_config.get_int_array("outlet_vessels")) {
      connections.push_back({junction_name, vessel_id_map[vessel_id]});
    }
    DEBUG_MSG("Created junction " << junction_name);
  }

  // Create Connections
  for (auto &connection : connections) {
    for (auto &ele1 : model->blocks) {
      for (auto &ele2 : model->blocks) {
        if ((ele1->name == std::get<0>(connection)) &&
            (ele2->name == std::get<1>(connection))) {
          model->nodes.push_back(
              new MODEL::Node(ele1->name + ":" + ele2->name));
          DEBUG_MSG("Created node " << model->nodes.back()->name);
          ele1->outlet_nodes.push_back(model->nodes.back());
          ele2->inlet_nodes.push_back(model->nodes.back());
          model->nodes.back()->setup_dofs(model->dofhandler);
        }
      }
    }
  }

  for (auto block : model->blocks) {
    if (block->inlet_nodes.size() == 0) {
      model->nodes.push_back(new MODEL::Node("inlet:" + block->name));
      block->inlet_nodes.push_back(model->nodes.back());
      model->nodes.back()->setup_dofs(model->dofhandler);
      DEBUG_MSG("Created node " << model->nodes.back()->name);
    }
    if (block->outlet_nodes.size() == 0) {
      model->nodes.push_back(new MODEL::Node(block->name + ":outlet"));
      block->outlet_nodes.push_back(model->nodes.back());
      model->nodes.back()->setup_dofs(model->dofhandler);
      DEBUG_MSG("Created node " << model->nodes.back()->name);
    }
  }
  // Setup degrees of freedom of the system
  for (auto &block : model->blocks) {
    block->setup_dofs(model->dofhandler);
  }

  for (size_t i = 0; i < model->dofhandler.size(); i++) {
    std::string var_name = model->dofhandler.variables[i];
    DEBUG_MSG("Variable " << i << ": " << var_name);
  }

  DEBUG_MSG("Number of parameters " << param_counter);

  DEBUG_MSG("Reading observations");
  std::vector<std::vector<T>> y_all;
  std::vector<std::vector<T>> dy_all;
  auto y_values = handler["y"];
  auto dy_values = handler["dy"];
  for (size_t i = 0; i < model->dofhandler.size(); i++) {
    std::string var_name = model->dofhandler.variables[i];
    DEBUG_MSG("Reading y values for variable " << var_name);
    y_all.push_back(y_values.get_double_array(var_name));
    DEBUG_MSG("Reading dy values for variable " << var_name);
    dy_all.push_back(dy_values.get_double_array(var_name));
  }

  int num_observations = y_all[0].size();
  DEBUG_MSG("Number of observations: " << num_observations);

  Eigen::SparseMatrix<T> X = Eigen::SparseMatrix<T>(
      num_observations * model->dofhandler.size(), param_counter);
  Eigen::Matrix<T, Eigen::Dynamic, 1> Y =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(num_observations *
                                                model->dofhandler.size());

  int num_dofs = model->dofhandler.size();
  for (size_t i = 0; i < num_observations; i++) {
    DEBUG_MSG("Assembling observation " << i);
    Eigen::Matrix<T, Eigen::Dynamic, 1> y =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(model->dofhandler.size());
    Eigen::Matrix<T, Eigen::Dynamic, 1> dy =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(model->dofhandler.size());
    for (size_t k = 0; k < model->dofhandler.size(); k++) {
      DEBUG_MSG("Loading DoF: " << k);
      y(k) = y_all[k][i];
      dy(k) = dy_all[k][i];
    }

    for (auto block : model->blocks) {
      DEBUG_MSG("Updating gradient of block " << block->id);
      block->update_gradient(X, Y, y, dy);

      for (size_t l = 0; l < block->global_eqn_ids.size(); l++) {
        block->global_eqn_ids[l] += num_dofs;
      }
    }
  }

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat = X.transpose() * X;
  Eigen::Matrix<T, Eigen::Dynamic, 1> alpha = mat.inverse() * X.transpose() * Y;

//   for (size_t i = 0; i < param_counter; i++) {
//     std::cout << alpha[i] << std::endl;
//   }
  for (auto block : model->blocks)
  {
    std::cout << block->name << std::endl;
    for (size_t j = 0; j < block->global_param_ids.size(); j++)
    {
        std::cout << alpha[block->global_param_ids[j]] << std::endl;
    }
    std::cout << "\n" << std::endl;
  }
}
