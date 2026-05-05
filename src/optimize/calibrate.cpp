// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "calibrate.h"

#include <set>

#include "LevenbergMarquardtOptimizer.h"
#include "SimulationParameters.h"

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
  bool calibrate_stenosis =
      calibration_parameters.value("calibrate_stenosis_coefficient", true);
  bool zero_capacitance =
      calibration_parameters.value("set_capacitance_to_zero", false);
  double lambda0 = calibration_parameters.value("initial_damping_factor", 1.0);

  // Optional list of parameter names to calibrate. Any parameter not listed is
  // held constant at its initial value from the input file. If the field is
  // absent or empty, all parameters are calibrated (legacy behavior).
  std::set<std::string> calibrate_names;
  bool calibrate_all = true;
  if (calibration_parameters.contains("calibrate")) {
    auto names =
        calibration_parameters["calibrate"].get<std::vector<std::string>>();
    if (!names.empty()) {
      calibrate_names.insert(names.begin(), names.end());
      calibrate_all = false;
    }
  }
  auto is_active = [&](const std::string& name) {
    return calibrate_all || calibrate_names.count(name) > 0;
  };

  int num_params = 3;
  if (calibrate_stenosis) {
    num_params = 4;
  }
  // Parameter names ordered to match BloodVessel::ParamId.
  const std::vector<std::string> bv_param_names = {
      "R_poiseuille", "C", "L", "stenosis_coefficient"};
  // Active parameter ids in alpha (ids of params actually optimized).
  std::vector<int> active_param_ids;

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

    // Create parameter IDs
    std::vector<int> param_ids;
    for (size_t k = 0; k < num_params; k++)
      param_ids.push_back(param_counter++);
    std::string block_type =
        vessel_config["zero_d_element_type"].get<std::string>();
    model.add_block(block_type, param_ids, vessel_name);
    vessel_id_map.insert({vessel_config["vessel_id"], vessel_name});
    DEBUG_MSG("Created vessel " << vessel_name);

    // Mark which of this block's parameters are active.
    for (size_t k = 0; k < num_params; k++) {
      if (is_active(bv_param_names[k])) {
        active_param_ids.push_back(param_ids[k]);
      }
    }

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
      std::vector<int> param_ids;
      for (size_t i = 0; i < (num_outlets * (num_params - 1)); i++)
        param_ids.push_back(param_counter++);
      model.add_block("BloodVesselJunction", param_ids, junction_name);

      // Mark which of this junction's per-outlet parameters are active.
      // Layout: [R0..Rn-1, L0..Ln-1, (S0..Sn-1)?]
      for (size_t i = 0; i < num_outlets; i++) {
        if (is_active("R_poiseuille"))
          active_param_ids.push_back(param_ids[i]);
      }
      for (size_t i = 0; i < num_outlets; i++) {
        if (is_active("L"))
          active_param_ids.push_back(param_ids[num_outlets + i]);
      }
      if (num_params > 3) {
        for (size_t i = 0; i < num_outlets; i++) {
          if (is_active("stenosis_coefficient"))
            active_param_ids.push_back(param_ids[2 * num_outlets + i]);
        }
      }
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
  DEBUG_MSG("Reading initial alpha");
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    DEBUG_MSG("Reading initial alpha for " << vessel_name);
    auto block = model.get_block(vessel_name);
    alpha[block->global_param_ids[0]] =
        vessel_config["zero_d_element_values"].value("R_poiseuille", 0.0);
    alpha[block->global_param_ids[1]] =
        vessel_config["zero_d_element_values"].value("C", 0.0);
    alpha[block->global_param_ids[2]] =
        vessel_config["zero_d_element_values"].value("L", 0.0);
    if (num_params > 3) {
      alpha[block->global_param_ids[3]] =
          vessel_config["zero_d_element_values"].value("stenosis_coefficient",
                                                       0.0);
    }
  }
  for (auto& junction_config : output_config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    DEBUG_MSG("Reading initial alpha for " << junction_name);
    auto block = model.get_block(junction_name);
    int num_outlets = block->outlet_nodes.size();

    if (num_outlets < 2) {
      continue;
    }

    for (size_t i = 0; i < num_outlets; i++) {
      alpha[block->global_param_ids[i]] = 0.0;
      alpha[block->global_param_ids[i + num_outlets]] = 0.0;
      if (num_params > 3) {
        alpha[block->global_param_ids[i + 2 * num_outlets]] = 0.0;
      }
    }
    if (junction_config["junction_type"] == "BloodVesselJunction") {
      auto resistance = junction_config["junction_values"]["R_poiseuille"]
                            .get<std::vector<double>>();
      auto inductance =
          junction_config["junction_values"]["L"].get<std::vector<double>>();
      auto stenosis_coeff =
          junction_config["junction_values"]["stenosis_coefficient"]
              .get<std::vector<double>>();
      for (size_t i = 0; i < num_outlets; i++) {
        alpha[block->global_param_ids[i]] = resistance[i];
        alpha[block->global_param_ids[i + num_outlets]] = inductance[i];
        if (num_params > 3) {
          alpha[block->global_param_ids[i + 2 * num_outlets]] =
              stenosis_coeff[i];
        }
      }
    }
  }

  // Run optimization
  DEBUG_MSG("Start optimization");
  DEBUG_MSG("Number of active parameters " << active_param_ids.size());
  if (active_param_ids.empty()) {
    throw std::runtime_error(
        "[svzerodcalibrator] No parameters selected for calibration. Either "
        "omit 'calibrate' from calibration_parameters or list at least one "
        "parameter name.");
  }
  auto lm_alg = LevenbergMarquardtOptimizer(
      &model, num_obs, param_counter, active_param_ids, lambda0, gradient_tol,
      increment_tol, max_iter);

  alpha = lm_alg.run(alpha, y_all, dy_all);

  // Write optimized simulation config file
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    auto block = model.get_block(vessel_name);
    double stenosis_coeff = 0.0;
    if (num_params > 3) {
      stenosis_coeff = alpha[block->global_param_ids[3]];
    }
    double c_value = 0.0;
    if (!zero_capacitance) {
      c_value = alpha[block->global_param_ids[1]];
    }
    vessel_config["zero_d_element_values"] = {
        {"R_poiseuille", alpha[block->global_param_ids[0]]},
        {"C", std::max(c_value, 0.0)},
        {"L", std::max(alpha[block->global_param_ids[2]], 0.0)},
        {"stenosis_coefficient", stenosis_coeff}};
  }
  for (auto& junction_config : output_config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    auto block = model.get_block(junction_name);
    int num_outlets = block->outlet_nodes.size();

    if (num_outlets < 2) {
      continue;
    }

    std::vector<double> r_values;
    for (size_t i = 0; i < num_outlets; i++) {
      r_values.push_back(alpha[block->global_param_ids[i]]);
    }
    std::vector<double> l_values;
    for (size_t i = 0; i < num_outlets; i++) {
      l_values.push_back(
          std::max(alpha[block->global_param_ids[i + num_outlets]], 0.0));
    }

    std::vector<double> ste_values;

    if (num_params > 3) {
      for (size_t i = 0; i < num_outlets; i++) {
        ste_values.push_back(
            alpha[block->global_param_ids[i + 2 * num_outlets]]);
      }
    } else {
      for (size_t i = 0; i < num_outlets; i++) {
        ste_values.push_back(0.0);
      }
    }

    junction_config["junction_type"] = "BloodVesselJunction";
    junction_config["junction_values"] = {{"R_poiseuille", r_values},
                                          {"L", l_values},
                                          {"stenosis_coefficient", ste_values}};
  }

  output_config.erase("y");
  output_config.erase("dy");
  output_config.erase("calibration_parameters");

  return output_config;
}
