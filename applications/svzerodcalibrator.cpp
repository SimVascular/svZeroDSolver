

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>

#include "helpers/debug.hpp"
#include "helpers/endswith.hpp"
#include "io/jsonhandler.hpp"
#include "model/model.hpp"
#include "optimize/levenbergmarquardtoptimizer.hpp"

#include <nlohmann/json.hpp>

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

  std::ifstream input_file_stream2(input_file);
  nlohmann::json output_config = nlohmann::json::parse(input_file_stream2);

  // Load model and configuration
  DEBUG_MSG("Read configuration");
  auto handler = IO::JsonHandler(config);

  auto model = std::shared_ptr<MODEL::Model<T>>(new MODEL::Model<T>());
  std::vector<std::tuple<std::string, std::string>> connections;
  std::vector<std::tuple<std::string, std::string>> inlet_connections;
  std::vector<std::tuple<std::string, std::string>> outlet_connections;

  auto calibration_parameters = handler["calibration_parameters"];
  bool calibrate_stenosis = calibration_parameters.get_bool("calibrate_stenosis_coefficient");
  bool zero_capacitance = calibration_parameters.get_bool("set_capacitance_to_zero");

  int num_params = 3;
  if (calibrate_stenosis){
    num_params = 4;
  }

  // Create vessels
  DEBUG_MSG("Load vessels");
  std::map<std::int64_t, std::string> vessel_id_map;
  auto vessels = handler["vessels"];
  int param_counter = 0;
  for (size_t i = 0; i < vessels.length(); i++) {
    auto vessel_config = vessels[i];
    std::string vessel_name = vessel_config.get_string("vessel_name");
    vessel_id_map.insert({vessel_config.get_int("vessel_id"), vessel_name});
    std::vector<int> param_ids;
    for (size_t k = 0; k < num_params; k++)
    {
        param_ids.push_back(param_counter++);
    }
    model->add_block(MODEL::BlockType::BLOODVESSEL,
                     param_ids,
                     vessel_name);
    DEBUG_MSG("Created vessel " << vessel_name);

    // Read connected boundary conditions
    if (vessel_config.has_key("boundary_conditions")) {
      auto vessel_bc_config = vessel_config["boundary_conditions"];
      if (vessel_bc_config.has_key("inlet")) {
        inlet_connections.push_back(
            {vessel_bc_config.get_string("inlet"), vessel_name});
      }
      if (vessel_bc_config.has_key("outlet")) {
        outlet_connections.push_back(
            {vessel_name, vessel_bc_config.get_string("outlet")});
      }
    }
  }

  // Create junctions
  auto junctions = handler["junctions"];
  for (size_t i = 0; i < junctions.length(); i++) {
    auto junction_config = junctions[i];
    auto junction_name = junction_config.get_string("junction_name");
    auto outlet_vessels = junction_config.get_int_array("outlet_vessels");
    int num_outlets = outlet_vessels.size();
    if (num_outlets == 1) {
      model->add_block(MODEL::BlockType::JUNCTION, {}, junction_name);
    } else {
      std::vector<int> param_ids;
      for (size_t i = 0; i < (num_outlets * num_params); i++) {
        param_ids.push_back(param_counter++);
      }
      model->add_block(MODEL::BlockType::BLOODVESSELJUNCTION, param_ids,
                       junction_name);
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
  DEBUG_MSG("Created connection");
  for (auto &connection : connections) {
    auto ele1 = model->get_block(std::get<0>(connection));
    auto ele2 = model->get_block(std::get<1>(connection));
    model->add_node({ele1}, {ele2}, ele1->get_name() + ":" + ele2->get_name());
  }

  for (auto &connection : inlet_connections) {
    auto ele = model->get_block(std::get<1>(connection));
    model->add_node({}, {ele}, static_cast<std::string>(std::get<0>(connection)) + ":" + ele->get_name());
  }

  for (auto &connection : outlet_connections) {
    auto ele = model->get_block(std::get<0>(connection));
    model->add_node({ele}, {}, ele->get_name() + ":" + static_cast<std::string>(std::get<1>(connection)));
  }

  model->finalize();

  for (size_t i = 0; i < model->dofhandler.size(); i++) {
    std::string var_name = model->dofhandler.variables[i];
    DEBUG_MSG("Variable " << i << ": " << var_name);
  }

  DEBUG_MSG("Number of parameters " << param_counter);
  int num_obs = 0;

  DEBUG_MSG("Reading observations");
  std::vector<std::vector<T>> y_all;
  std::vector<std::vector<T>> dy_all;
  auto y_values = handler["y"];
  auto dy_values = handler["dy"];
  for (size_t i = 0; i < model->dofhandler.get_num_variables(); i++) {

    std::string var_name = model->dofhandler.variables[i];
    DEBUG_MSG("Reading observations for variable " << var_name);
    auto y_array = y_values.get_double_array(var_name);
    auto dy_array = dy_values.get_double_array(var_name);
    num_obs = y_array.size();
    if (i==0)
    {
      y_all.resize(num_obs);
      dy_all.resize(num_obs);
    }
    for (size_t j = 0; j < num_obs; j++)
    {
      y_all[j].push_back(y_array[j]);
      dy_all[j].push_back(dy_array[j]);
    }
  }
  DEBUG_MSG("Number of observations: " << num_obs);

  // =====================================
  // Setup alpha
  // =====================================
  Eigen::Matrix<T, Eigen::Dynamic, 1> alpha =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(param_counter);
  DEBUG_MSG("Reading initial alpha");
  for (auto &vessel_config : output_config["vessels"])
  {
    std::string vessel_name = vessel_config["vessel_name"];
    DEBUG_MSG("Reading initial alpha for " << vessel_name);
    auto block = model->get_block(vessel_name);
    alpha[block->global_param_ids[0]] = vessel_config["zero_d_element_values"].value("R_poiseuille", 0.0);
    alpha[block->global_param_ids[1]] = vessel_config["zero_d_element_values"].value("C", 0.0);
    alpha[block->global_param_ids[2]] = vessel_config["zero_d_element_values"].value("L", 0.0);
    if (num_params > 3)
    {
      alpha[block->global_param_ids[3]] = vessel_config["zero_d_element_values"].value("stenosis_coefficient", 0.0);
    }
  }
  for (auto &junction_config : output_config["junctions"])
  {

    std::string junction_name = junction_config["junction_name"];
    DEBUG_MSG("Reading initial alpha for " << junction_name);
    auto block = model->get_block(junction_name);
    int num_outlets = block->outlet_nodes.size();

    if (num_outlets < 2)
    {
        continue;
    }

    for (size_t i = 0; i < num_outlets; i++)
    {
        alpha[block->global_param_ids[i]] = 0.0;
        alpha[block->global_param_ids[i+num_outlets]] = 0.0;
        alpha[block->global_param_ids[i+2*num_outlets]] = 0.0;
        if (num_params > 3) {
          alpha[block->global_param_ids[i+3*num_outlets]] = 0.0;
        }
    }
    if (junction_config["junction_type"] == "BloodVesselJunction")
    {
      auto resistance = junction_config["junction_values"]["R_poiseuille"].get<std::vector<double>>();
      auto capacitance = junction_config["junction_values"]["C"].get<std::vector<double>>();
      auto inductance = junction_config["junction_values"]["L"].get<std::vector<double>>();
      auto stenosis_coeff = junction_config["junction_values"]["stenosis_coefficient"].get<std::vector<double>>();
      for (size_t i = 0; i < num_outlets; i++)
      {
          alpha[block->global_param_ids[i]] = resistance[i];
          alpha[block->global_param_ids[i+num_outlets]] = capacitance[i];
          alpha[block->global_param_ids[i+2*num_outlets]] = inductance[i];
          if (num_params > 3) {
            alpha[block->global_param_ids[i+3*num_outlets]] = stenosis_coeff[i];
          }
      }
    }
  }


  // =====================================
  // Levenberg Marquardt
  // =====================================
  DEBUG_MSG("Start optimization");
  auto lm_alg = OPT::LevenbergMarquardtOptimizer(model.get(), num_obs, param_counter, 1.0);
  
  alpha = lm_alg.run(alpha, y_all, dy_all);

  // =====================================
  // Write output
  // =====================================
  for (auto &vessel_config : output_config["vessels"])
  {
    std::string vessel_name = vessel_config["vessel_name"];
    auto block = model->get_block(vessel_name);
    T stenosis_coeff = 0.0;
    if (num_params > 3)
    {
      stenosis_coeff = alpha[block->global_param_ids[3]];
    }
    T c_value = 0.0;
    if (!zero_capacitance)
    {
      c_value = alpha[block->global_param_ids[1]];
    }
    vessel_config["zero_d_element_values"] = {
        {"R_poiseuille", alpha[block->global_param_ids[0]]},
        {"C", std::max(c_value, 0.0)},
        {"L", std::max(alpha[block->global_param_ids[2]], 0.0)},
        {"stenosis_coefficient", stenosis_coeff}
    };
  }

  for (auto &junction_config : output_config["junctions"])
  {

    std::string junction_name = junction_config["junction_name"];
    auto block = model->get_block(junction_name);
    int num_outlets = block->outlet_nodes.size();

    if (num_outlets < 2)
    {
        continue;
    }

    std::vector<T> r_values;
    for (size_t i = 0; i < num_outlets; i++)
    {
        r_values.push_back(alpha[block->global_param_ids[i]]);
    }
    std::vector<T> c_values;
    if (zero_capacitance)
    {
      for (size_t i = 0; i < num_outlets; i++)
      {
          c_values.push_back(0.0);
      }
    }
    else {
      for (size_t i = 0; i < num_outlets; i++)
      {
          c_values.push_back(std::max(alpha[block->global_param_ids[i+num_outlets]], 0.0));
      }
    }
    std::vector<T> l_values;
    for (size_t i = 0; i < num_outlets; i++)
    {
        l_values.push_back(std::max(alpha[block->global_param_ids[i+2*num_outlets]],0.0));
    }
    std::vector<T> ste_values;
    if (num_params > 3)
    {
        for (size_t i = 0; i < num_outlets; i++){
          ste_values.push_back(alpha[block->global_param_ids[i+3*num_outlets]]);
        }
    }
    else {
      for (size_t i = 0; i < num_outlets; i++)
      {
          ste_values.push_back(0.0);
      }
    }
    junction_config["junction_type"] = "BloodVesselJunction";
    junction_config["junction_values"] = {
        {"R_poiseuille", r_values},
        {"C", c_values},
        {"L", l_values},
        {"stenosis_coefficient", ste_values}
    };
  }

  output_config.erase("y");
  output_config.erase("dy");
  output_config.erase("calibration_parameters");

  std::ofstream o(output_file);
  o << std::setw(4) << output_config << std::endl;
}
