// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/**
 * @file configreader.hpp
 * @brief IO::ConfigReader source file
 */
#ifndef SVZERODSOLVER_IO_CONFIGREADER_HPP_
#define SVZERODSOLVER_IO_CONFIGREADER_HPP_

#include <list>
#include <stdexcept>
#include <string>

#include "../helpers/debug.hpp"
#include "../helpers/startswith.hpp"
#include "../model/bloodvessel.hpp"
#include "../model/flowreferencebc.hpp"
#include "../model/junction.hpp"
#include "../model/model.hpp"
#include "../model/node.hpp"
#include "../model/openloopcoronarybc.hpp"
#include "../model/pressurereferencebc.hpp"
#include "../model/resistancebc.hpp"
#include "../model/resistivejunction.hpp"
#include "../model/windkesselbc.hpp"
#include "simdjson.h"

namespace IO {

/**
 * @brief svZeroDSolver configuration reader.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class ConfigReader {
 public:
  /**
   * @brief Construct a new Config Reader object
   */
  ConfigReader();

  /**
   * @brief Destroy the Config Reader object
   */
  ~ConfigReader();

  /**
   * @brief Create model from configuration file
   *
   * @return Model
   */
  void load(std::string &specifier);

  MODEL::Model<T> model;            ///< Simulation model
  ALGEBRA::State<T> initial_state;  ///< Initial state

  T sim_cardiac_cycle_period;  ///< Cardiac cycle period
  T sim_time_step_size;        ///< Simulation time step size
  T sim_abs_tol;               ///< Absolute tolerance for simulation

  int sim_num_cycles;      ///< Number of cardiac cycles to simulate
  int sim_pts_per_cycle;   ///< Number of time steps per cardiac cycle
  int sim_num_time_steps;  ///< Total number of time steps
  int sim_nliter;          ///< Maximum number of non-linear iterations in time
                           ///< integration
  int output_interval;     ///< Interval of writing output

  bool sim_steady_initial;  ///< Start from steady solution
  bool
      output_variable_based;  ///< Output variable based instead of vessel based
  bool output_mean_only;      ///< Output only the mean value
  bool output_derivative;     ///< Output derivatives
  bool output_last_cycle_only;  ///< Output only the last cardiac cycle
};

template <typename T>
ConfigReader<T>::ConfigReader() {}

template <typename T>
ConfigReader<T>::~ConfigReader() {}

template <typename T>
void ConfigReader<T>::load(std::string &specifier) {
  // Create iterator for json configuration
  simdjson::ondemand::parser parser;
  simdjson::padded_string string;
  if (HELPERS::startswith(specifier, "{")) {
    string = simdjson::padded_string(specifier);
  } else {
    string = simdjson::padded_string::load(specifier);
  }
  auto config = parser.iterate(string);

  // Load simulation paramaters
  auto sim_params = config["simulation_parameters"];
  sim_num_cycles = sim_params["number_of_cardiac_cycles"].get_int64();
  sim_pts_per_cycle =
      sim_params["number_of_time_pts_per_cardiac_cycle"].get_int64();
  sim_num_time_steps = (sim_pts_per_cycle - 1.0) * sim_num_cycles + 1.0;
  try {
    sim_abs_tol = sim_params["absolute_tolerance"].get_double();
  } catch (simdjson::simdjson_error) {
    sim_abs_tol = 1e-8;
  }
  try {
    sim_nliter = sim_params["maximum_nonlinear_iterations"].get_int64();
  } catch (simdjson::simdjson_error) {
    sim_nliter = 30;
  }
  try {
    sim_steady_initial = sim_params["steady_initial"].get_bool();
  } catch (simdjson::simdjson_error) {
    sim_steady_initial = true;
  }
  try {
    output_variable_based = sim_params["output_variable_based"].get_bool();
  } catch (simdjson::simdjson_error) {
    output_variable_based = false;
  }
  try {
    output_interval = sim_params["output_interval"].get_int64();
  } catch (simdjson::simdjson_error) {
    output_interval = 1;
  }
  try {
    output_mean_only = sim_params["output_mean_only"].get_bool();
  } catch (simdjson::simdjson_error) {
    output_mean_only = false;
  }
  try {
    output_derivative = sim_params["output_derivative"].get_bool();
  } catch (simdjson::simdjson_error) {
    output_derivative = false;
  }
  try {
    output_last_cycle_only = sim_params["output_last_cycle_only"].get_bool();
  } catch (simdjson::simdjson_error) {
    output_last_cycle_only = false;
  }

  // Create list to store block connections while generating blocks
  std::vector<std::tuple<std::string_view, std::string_view>> connections;

  // Create vessels
  std::map<std::int64_t, std::string_view> vessel_id_map;
  for (auto vessel_config : config["vessels"]) {
    auto vessel_values = vessel_config["zero_d_element_values"];
    std::string_view vessel_name = vessel_config["vessel_name"].get_string();
    vessel_id_map.insert({vessel_config["vessel_id"].get_int64(), vessel_name});
    if (std::string_view(vessel_config["zero_d_element_type"].get_string()) ==
        "BloodVessel") {
      T R = vessel_values["R_poiseuille"].get_double();
      T C;
      try {
        C = vessel_values["C"].get_double();
      } catch (simdjson::simdjson_error) {
        C = 0.0;
      }
      T L;
      try {
        L = vessel_values["L"].get_double();
      } catch (simdjson::simdjson_error) {
        L = 0.0;
      }
      T stenosis_coefficient;
      try {
        stenosis_coefficient =
            vessel_values["stenosis_coefficient"].get_double();
      } catch (simdjson::simdjson_error) {
        stenosis_coefficient = 0.0;
      }
      model.blocks.push_back(new MODEL::BloodVessel<T>(
          R = R, C = C, L = L, stenosis_coefficient = stenosis_coefficient,
          static_cast<std::string>(vessel_name)));
      DEBUG_MSG("Created vessel " << vessel_name);
    } else {
      throw std::invalid_argument("Unknown vessel type");
    }

    // Read connected boundary condtitions
    try {
      auto inlet_bc = vessel_config["boundary_conditions"]["inlet"];
      connections.push_back({inlet_bc.get_string(), vessel_name});
    } catch (simdjson::simdjson_error) {
    }

    try {
      auto outlet_bc = vessel_config["boundary_conditions"]["outlet"];
      connections.push_back({vessel_name, outlet_bc.get_string()});
    } catch (simdjson::simdjson_error) {
    }
  }

  // Create boundary conditions
  for (auto bc_config : config["boundary_conditions"]) {
    std::string_view bc_type = bc_config["bc_type"];
    std::string_view bc_name = bc_config["bc_name"];
    auto bc_values = bc_config["bc_values"];

    std::vector<T> t;
    try {
      for (auto x : bc_values["t"]) {
        t.push_back(x.get_double());
      }
    } catch (simdjson::simdjson_error) {
      t.push_back(0.0);
    };

    if (bc_type == "RCR") {
      T Rp = bc_values["Rp"], C = bc_values["C"], Rd = bc_values["Rd"],
        Pd = bc_values["Pd"];
      model.blocks.push_back(new MODEL::WindkesselBC<T>(
          Rp = Rp, C = C, Rd = Rd, Pd = Pd, static_cast<std::string>(bc_name)));
      DEBUG_MSG("Created boundary condition " << bc_name);
    } else if (bc_type == "FLOW") {
      std::vector<T> Q;
      try {
        for (auto x : bc_values["Q"].get_array()) {
          Q.push_back(x.get_double());
        }
      } catch (simdjson::simdjson_error) {
        Q.push_back(bc_values["Q"].get_double());
      }

      MODEL::TimeDependentParameter q_param(t, Q);
      if (q_param.isconstant == false) {
        sim_cardiac_cycle_period = q_param.cycle_period;
      }
      model.blocks.push_back(new MODEL::FlowReferenceBC<T>(
          q_param, static_cast<std::string>(bc_name)));
      DEBUG_MSG("Created boundary condition " << bc_name);

    } else if (bc_type == "RESISTANCE") {
      std::vector<T> R;
      try {
        for (auto x : bc_values["R"].get_array()) {
          R.push_back(x.get_double());
        }
      } catch (simdjson::simdjson_error) {
        R.push_back(bc_values["R"].get_double());
      }
      MODEL::TimeDependentParameter r_param(t, R);

      std::vector<T> Pd;
      try {
        for (auto x : bc_values["Pd"].get_array()) {
          Pd.push_back(x.get_double());
        }
      } catch (simdjson::simdjson_error) {
        Pd.push_back(bc_values["Pd"].get_double());
      }
      MODEL::TimeDependentParameter pd_param(t, Pd);
      model.blocks.push_back(new MODEL::ResistanceBC<T>(
          r_param, pd_param, static_cast<std::string>(bc_name)));
      DEBUG_MSG("Created boundary condition " << bc_name);
    } else if (bc_type == "PRESSURE") {
      std::vector<T> P;
      try {
        for (auto x : bc_values["P"].get_array()) {
          P.push_back(x.get_double());
        }
      } catch (simdjson::simdjson_error) {
        P.push_back(bc_values["P"].get_double());
      }
      MODEL::TimeDependentParameter p_param(t, P);
      if (p_param.isconstant == false) {
        sim_cardiac_cycle_period = p_param.cycle_period;
      }
      model.blocks.push_back(new MODEL::PressureReferenceBC<T>(
          p_param, static_cast<std::string>(bc_name)));
      DEBUG_MSG("Created boundary condition " << bc_name);
    } else if (bc_type == "CORONARY") {
      T Ra = bc_values["Ra1"], Ram = bc_values["Ra2"], Rv = bc_values["Rv1"],
        Ca = bc_values["Ca"], Cim = bc_values["Cc"];
      auto Pim_json = bc_values["Pim"], Pv_json = bc_values["P_v"];

      std::vector<T> Pim;
      try {
        for (auto x : bc_values["Pim"].get_array()) {
          Pim.push_back(x.get_double());
        }
      } catch (simdjson::simdjson_error) {
        Pim.push_back(bc_values["Pim"].get_double());
      }
      MODEL::TimeDependentParameter pim_param(t, Pim);
      std::vector<T> P_v;
      try {
        for (auto x : bc_values["P_v"].get_array()) {
          P_v.push_back(x.get_double());
        }
      } catch (simdjson::simdjson_error) {
        P_v.push_back(bc_values["P_v"].get_double());
      }
      MODEL::TimeDependentParameter pv_param(t, P_v);

      model.blocks.push_back(new MODEL::OpenLoopCoronaryBC<T>(
          Ra = Ra, Ram = Ram, Rv = Rv, Ca = Ca, Cim = Cim, pim_param, pv_param,
          static_cast<std::string>(bc_name)));
      DEBUG_MSG("Created boundary condition " << bc_name);
    } else {
      throw std::invalid_argument("Unknown boundary condition type");
    }
  }

  // Create junctions
  for (auto junction_config : config["junctions"]) {
    std::string_view j_type = junction_config["junction_type"].get_string();
    std::string_view junction_name =
        junction_config["junction_name"].get_string();
    if ((j_type == "NORMAL_JUNCTION") || (j_type == "internal_junction")) {
      model.blocks.push_back(
          new MODEL::Junction<T>(static_cast<std::string>(junction_name)));
    } else if (j_type == "resistive_junction") {
      std::vector<T> R;
      for (auto x : junction_config["junction_values"]["R"].get_array()) {
        R.push_back(x.get_double());
      }
      model.blocks.push_back(new MODEL::ResistiveJunction<T>(
          R, static_cast<std::string>(junction_name)));
    } else {
      throw std::invalid_argument("Unknown junction type");
    }
    // Check for connections to inlet and outlet vessels and append to
    // connections list
    for (auto inlet_vessel : junction_config["inlet_vessels"]) {
      connections.push_back(
          {vessel_id_map[inlet_vessel.get_int64()], junction_name});
    }
    for (auto outlet_vessel : junction_config["outlet_vessels"]) {
      connections.push_back(
          {junction_name, vessel_id_map[outlet_vessel.get_int64()]});
    }
  }

  // Create Connections
  for (auto &connection : connections) {
    for (auto &ele1 : model.blocks) {
      for (auto &ele2 : model.blocks) {
        if ((ele1->name == std::get<0>(connection)) &&
            (ele2->name == std::get<1>(connection))) {
          model.nodes.push_back(new MODEL::Node(ele1->name + ":" + ele2->name));
          DEBUG_MSG("Created node " << model.nodes.back()->name);
          ele1->outlet_nodes.push_back(model.nodes.back());
          ele2->inlet_nodes.push_back(model.nodes.back());
          model.nodes.back()->setup_dofs(model.dofhandler);
        }
      }
    }
  }

  // Setup degrees of freedom of the system
  for (auto &block : model.blocks) {
    block->setup_dofs(model.dofhandler);
  }

  // Read initial condition
  initial_state = ALGEBRA::State<T>::Zero(model.dofhandler.size());
  try {
    auto initial_condition = config["initial_condition"].value();
    for (size_t i = 0; i < model.dofhandler.size(); i++) {
      try {
        initial_state.y[i] = initial_condition[model.dofhandler.variables[i]];
      } catch (simdjson::simdjson_error) {
        std::cout << "Initial condition found, but no value for "
                  << model.dofhandler.variables[i] << std::endl;
      }
    }
  } catch (simdjson::simdjson_error) {
  }

  // Calculate time step size
  sim_time_step_size = sim_cardiac_cycle_period / (T(sim_pts_per_cycle) - 1.0);
}

}  // namespace IO

#endif  // SVZERODSOLVER_IO_CONFIGREADER_HPP_