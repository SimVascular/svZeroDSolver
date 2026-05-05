// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "calibrate.h"

#include <set>

#include "LevenbergMarquardtOptimizer.h"
#include "SimulationParameters.h"

namespace {

// Number of alpha slots a block exposes to the calibrator. Scalar blocks
// (``input_params_list == false``) contribute one slot per named input
// parameter. List blocks contribute ``input_params.size() * stride`` slots,
// where ``stride`` is determined by the block's connectivity (e.g. number of
// outlets for BloodVesselJunction).
int num_param_slots(const Block& block, int stride) {
  int n = static_cast<int>(block.input_params.size());
  return block.input_params_list ? n * stride : n;
}

// Stride for a list-style block. The only list block currently in the model
// (``BloodVesselJunction``) groups its parameters per outlet, so the stride
// equals the number of outlets. For scalar blocks the stride is 1.
int param_stride(const Block& block, int num_outlets) {
  return block.input_params_list ? num_outlets : 1;
}

}  // namespace

nlohmann::json calibrate(const nlohmann::json& config) {
  auto output_config = nlohmann::json(config);

  // Read calibration parameters
  DEBUG_MSG("Parse calibration parameters");
  auto const& calibration_parameters = config["calibration_parameters"];
  double gradient_tol =
      calibration_parameters.value("tolerance_gradient", 1e-5);
  double increment_tol =
      calibration_parameters.value("tolerance_increment", 1e-10);
  int max_iter = calibration_parameters.value("maximum_iterations", 100);
  double lambda0 = calibration_parameters.value("initial_damping_factor", 1.0);

  // Resolve the set of parameter names to calibrate for a given block from
  // its own ``calibrate`` field. The field is mandatory: a block without a
  // ``calibrate`` field has no parameters optimized (every parameter of that
  // block is held at its input value). An explicit empty list is equivalent
  // to omitting the field.
  auto resolve_calibrate_set = [](const nlohmann::json& block_config) {
    std::set<std::string> names;
    if (!block_config.contains("calibrate")) return names;
    auto list = block_config["calibrate"].get<std::vector<std::string>>();
    names.insert(list.begin(), list.end());
    return names;
  };

  // Append the active alpha indices contributed by a single block to
  // ``active_param_ids``. Walks the block's ``input_params`` and includes
  // each parameter (or, for list blocks, each per-outlet copy of it) whose
  // name appears in the resolved calibrate filter.
  std::vector<int> active_param_ids;
  auto register_active = [&](const Block& block,
                             const std::vector<int>& param_ids,
                             const std::set<std::string>& set) {
    int stride = static_cast<int>(param_ids.size()) /
                 static_cast<int>(block.input_params.size());
    if (!block.input_params_list) stride = 1;
    for (size_t i = 0; i < block.input_params.size(); i++) {
      const std::string& name = block.input_params[i].first;
      if (set.count(name) == 0) continue;
      for (int s = 0; s < stride; s++) {
        active_param_ids.push_back(param_ids[i * stride + s]);
      }
    }
  };

  // Setup model
  auto model = Model();
  std::vector<std::tuple<std::string, std::string>> connections;
  std::vector<std::tuple<std::string, std::string>> inlet_connections;
  std::vector<std::tuple<std::string, std::string>> outlet_connections;

  // Create vessels
  DEBUG_MSG("Load vessels");
  std::map<std::int64_t, std::string> vessel_id_map;
  int param_counter = 0;
  for (auto const& vessel_config : config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    std::string block_type =
        vessel_config["zero_d_element_type"].get<std::string>();

    // Instantiate the block so we can introspect its parameter names.
    auto* block = model.create_block(block_type);
    int num_slots = num_param_slots(*block, /*stride=*/1);
    std::vector<int> param_ids;
    for (int k = 0; k < num_slots; k++) param_ids.push_back(param_counter++);
    model.add_block(block, vessel_name, param_ids);
    vessel_id_map.insert({vessel_config["vessel_id"], vessel_name});
    DEBUG_MSG("Created vessel " << vessel_name);

    register_active(*block, param_ids, resolve_calibrate_set(vessel_config));

    // Read connected boundary conditions
    if (vessel_config.contains("boundary_conditions")) {
      auto const& vessel_bc_config = vessel_config["boundary_conditions"];
      if (vessel_bc_config.contains("inlet")) {
        inlet_connections.push_back({vessel_bc_config["inlet"], vessel_name});
      }
      if (vessel_bc_config.contains("outlet")) {
        outlet_connections.push_back({vessel_name, vessel_bc_config["outlet"]});
      }
    }
  }

  // Create junctions
  for (auto const& junction_config : config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    auto const& outlet_vessels = junction_config["outlet_vessels"];
    int num_outlets = outlet_vessels.size();

    if (num_outlets == 1) {
      model.add_block("NORMAL_JUNCTION", {}, junction_name);
    } else {
      auto* block = model.create_block("BloodVesselJunction");
      int num_slots = num_param_slots(*block, num_outlets);
      std::vector<int> param_ids;
      for (int k = 0; k < num_slots; k++) param_ids.push_back(param_counter++);
      model.add_block(block, junction_name, param_ids);

      register_active(*block, param_ids,
                      resolve_calibrate_set(junction_config));
    }

    // Check for connections to inlet and outlet vessels and append to
    // connections list
    for (auto vessel_id : junction_config["inlet_vessels"]) {
      connections.push_back({vessel_id_map[vessel_id], junction_name});
    }

    for (auto vessel_id : outlet_vessels) {
      connections.push_back({junction_name, vessel_id_map[vessel_id]});
    }
    DEBUG_MSG("Created junction " << junction_name);
  }

  // Create Connections
  DEBUG_MSG("Created connection");
  for (auto& connection : connections) {
    auto ele1 = model.get_block(std::get<0>(connection));
    auto ele2 = model.get_block(std::get<1>(connection));
    model.add_node({ele1}, {ele2}, ele1->get_name() + ":" + ele2->get_name());
  }
  for (auto& connection : inlet_connections) {
    auto ele = model.get_block(std::get<1>(connection));
    model.add_node({}, {ele}, std::get<0>(connection) + ":" + ele->get_name());
  }
  for (auto& connection : outlet_connections) {
    auto ele = model.get_block(std::get<0>(connection));
    model.add_node({ele}, {}, ele->get_name() + ":" + std::get<1>(connection));
  }

  // Finalize model
  model.finalize();

  DEBUG_MSG("Number of parameters " << param_counter);

  // Read observations
  DEBUG_MSG("Reading observations");
  int num_obs = 0;
  std::vector<std::vector<double>> y_all;
  std::vector<std::vector<double>> dy_all;
  auto y_values = config["y"];
  auto dy_values = config["dy"];
  for (size_t i = 0; i < model.dofhandler.get_num_variables(); i++) {
    std::string var_name = model.dofhandler.variables[i];
    DEBUG_MSG("Reading observations for variable " << var_name);
    if (!y_values.contains(var_name)) {
      std::cout << "ERROR: Missing y observation for '" << var_name << "'"
                << std::endl;
      exit(1);
    }
    if (!dy_values.contains(var_name)) {
      std::cout << "ERROR: Missing dy observation for '" << var_name << "'"
                << std::endl;
      exit(1);
    }
    auto y_array = y_values[var_name].get<std::vector<double>>();
    auto dy_array = dy_values[var_name].get<std::vector<double>>();
    num_obs = y_array.size();
    if (i == 0) {
      y_all.resize(num_obs);
      dy_all.resize(num_obs);
    }
    for (size_t j = 0; j < num_obs; j++) {
      y_all[j].push_back(y_array[j]);
      dy_all[j].push_back(dy_array[j]);
    }
  }
  DEBUG_MSG("Number of observations: " << num_obs);

  // Setup start parameter vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(param_counter);

  // Initialize alpha from a JSON values object using the block's own
  // ``input_params`` to map names to slot indices. Missing names default to
  // 0 (already set by the zero-init above).
  auto init_alpha_for_block = [&](const Block& block,
                                  const nlohmann::json& values) {
    int total = static_cast<int>(block.global_param_ids.size());
    int stride = param_stride(
        block, total / static_cast<int>(block.input_params.size()));
    for (size_t i = 0; i < block.input_params.size(); i++) {
      const std::string& name = block.input_params[i].first;
      if (!values.contains(name)) continue;
      if (block.input_params_list) {
        auto arr = values[name].get<std::vector<double>>();
        for (int s = 0; s < stride && s < static_cast<int>(arr.size()); s++) {
          alpha[block.global_param_ids[i * stride + s]] = arr[s];
        }
      } else {
        alpha[block.global_param_ids[i]] = values[name].get<double>();
      }
    }
  };

  DEBUG_MSG("Reading initial alpha");
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    DEBUG_MSG("Reading initial alpha for " << vessel_name);
    auto block = model.get_block(vessel_name);
    if (vessel_config.contains("zero_d_element_values")) {
      init_alpha_for_block(*block, vessel_config["zero_d_element_values"]);
    }
  }
  for (auto& junction_config : output_config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    DEBUG_MSG("Reading initial alpha for " << junction_name);
    auto block = model.get_block(junction_name);
    if (block->global_param_ids.empty()) continue;
    if (junction_config.contains("junction_values")) {
      init_alpha_for_block(*block, junction_config["junction_values"]);
    }
  }

  // Run optimization
  DEBUG_MSG("Start optimization");
  DEBUG_MSG("Number of active parameters " << active_param_ids.size());
  if (active_param_ids.empty()) {
    throw std::runtime_error(
        "[svzerodcalibrator] No parameters selected for calibration. Add a "
        "'calibrate' field listing parameter names to at least one vessel or "
        "junction.");
  }
  auto lm_alg = LevenbergMarquardtOptimizer(
      &model, num_obs, param_counter, active_param_ids, lambda0, gradient_tol,
      increment_tol, max_iter);

  alpha = lm_alg.run(alpha, y_all, dy_all);

  // Build a JSON values object for a block by reading optimized alpha values
  // out using the block's own ``input_params``.
  auto write_alpha_for_block = [&](const Block& block) -> nlohmann::json {
    nlohmann::json values = nlohmann::json::object();
    int total = static_cast<int>(block.global_param_ids.size());
    int stride = param_stride(
        block, total / static_cast<int>(block.input_params.size()));
    for (size_t i = 0; i < block.input_params.size(); i++) {
      const std::string& name = block.input_params[i].first;
      if (block.input_params_list) {
        std::vector<double> arr;
        for (int s = 0; s < stride; s++) {
          arr.push_back(alpha[block.global_param_ids[i * stride + s]]);
        }
        values[name] = arr;
      } else {
        values[name] = alpha[block.global_param_ids[i]];
      }
    }
    return values;
  };

  // Write optimized simulation config file
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    auto block = model.get_block(vessel_name);
    vessel_config["zero_d_element_values"] = write_alpha_for_block(*block);
  }
  for (auto& junction_config : output_config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    auto block = model.get_block(junction_name);
    if (block->global_param_ids.empty()) continue;
    junction_config["junction_type"] = "BloodVesselJunction";
    junction_config["junction_values"] = write_alpha_for_block(*block);
  }

  output_config.erase("y");
  output_config.erase("dy");
  output_config.erase("calibration_parameters");

  return output_config;
}
