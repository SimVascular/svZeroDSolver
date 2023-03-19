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
#include "../model/block.hpp"
#include "../model/bloodvessel.hpp"
#include "../model/bloodvesseljunction.hpp"
#include "../model/closedloopRCRbc.hpp"
#include "../model/closedloopcoronarybc.hpp"
#include "../model/closedloopheartpulmonary.hpp"
#include "../model/flowreferencebc.hpp"
#include "../model/junction.hpp"
#include "../model/model.hpp"
#include "../model/node.hpp"
#include "../model/openloopcoronarybc.hpp"
#include "../model/pressurereferencebc.hpp"
#include "../model/resistancebc.hpp"
#include "../model/resistivejunction.hpp"
#include "../model/windkesselbc.hpp"
#include "./jsonhandler.hpp"

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
   * @brief Create model from configuration
   *
   * @param handler Configuration handler
   */
  void load(JsonHandler &handler);

  MODEL::Model<T> *model;           ///< Simulation model
  ALGEBRA::State<T> initial_state;  ///< Initial state

  T sim_cardiac_cycle_period;  ///< Cardiac cycle period.
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

  bool sim_coupled;  ///< Running 0D simulation coupled with external solver
  T sim_external_step_size;  ///< Step size of external solver if running
                             ///< coupled
};

template <typename T>
ConfigReader<T>::ConfigReader() {}

template <typename T>
ConfigReader<T>::~ConfigReader() {}

template <typename T>
void ConfigReader<T>::load(JsonHandler &handler) {
  DEBUG_MSG("Start loading configuration");
  // Initialize model pointer
  model = new MODEL::Model<T>();

  // Load simulation parameters
  DEBUG_MSG("Load simulation parameters");
  auto sim_params = handler["simulation_parameters"];
  sim_coupled = sim_params.get_bool("coupled_simulation", false);
  if (!sim_coupled) {
    sim_num_cycles = sim_params.get_int("number_of_cardiac_cycles");
    sim_pts_per_cycle =
        sim_params.get_int("number_of_time_pts_per_cardiac_cycle");
    sim_num_time_steps = (sim_pts_per_cycle - 1) * sim_num_cycles + 1;
    sim_external_step_size = 0.0;
  } else {
    sim_num_cycles = 1;
    sim_num_time_steps = sim_params.get_int("number_of_time_pts");
    sim_pts_per_cycle = sim_num_time_steps;
    sim_external_step_size = sim_params.get_double("external_step_size", 0.1);
  }
  sim_cardiac_cycle_period = -1.0;  // Negative value indicates this has not
                                    // been read from config file yet.
  sim_abs_tol = sim_params.get_double("absolute_tolerance", 1e-8);
  sim_nliter = sim_params.get_int("maximum_nonlinear_iterations", 30);
  sim_steady_initial = sim_params.get_bool("steady_initial", true);
  output_variable_based = sim_params.get_bool("output_variable_based", false);
  output_interval = sim_params.get_int("output_interval", 1);
  output_mean_only = sim_params.get_bool("output_mean_only", false);
  output_derivative = sim_params.get_bool("output_derivative", false);
  output_last_cycle_only = sim_params.get_bool("output_last_cycle_only", false);

  // Create list to store block connections while generating blocks
  std::vector<std::tuple<std::string_view, std::string_view>> connections;
  int block_count = 0;

  // Create vessels
  DEBUG_MSG("Load vessels");
  std::map<std::int64_t, std::string_view> vessel_id_map;
  auto vessels = handler["vessels"];
  for (size_t i = 0; i < vessels.length(); i++) {
    auto vessel_config = vessels[i];
    auto vessel_values = vessel_config["zero_d_element_values"];
    std::string_view vessel_name = vessel_config.get_string("vessel_name");
    vessel_id_map.insert({vessel_config.get_int("vessel_id"), vessel_name});
    if (vessel_config.get_string("zero_d_element_type") == "BloodVessel") {
      T R = vessel_values.get_double("R_poiseuille");
      T C = vessel_values.get_double("C", 0.0);
      T L = vessel_values.get_double("L", 0.0);
      T stenosis_coefficient =
          vessel_values.get_double("stenosis_coefficient", 0.0);
      model->blocks.push_back(new MODEL::BloodVessel<T>(
          R = R, C = C, L = L, stenosis_coefficient = stenosis_coefficient,
          static_cast<std::string>(vessel_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(vessel_name), block_count});
      block_count++;
      DEBUG_MSG("Created vessel " << vessel_name);
    } else {
      throw std::invalid_argument("Unknown vessel type");
    }

    // Read connected boundary conditions
    if (vessel_config.has_key("boundary_conditions")) {
      auto vessel_bc_config = vessel_config["boundary_conditions"];
      if (vessel_bc_config.has_key("inlet")) {
        connections.push_back(
            {vessel_bc_config.get_string("inlet"), vessel_name});
      }
      if (vessel_bc_config.has_key("outlet")) {
        connections.push_back(
            {vessel_name, vessel_bc_config.get_string("outlet")});
      }
    }
  }

  // Create map for boundary conditions to boundary condition type
  DEBUG_MSG("Create BC name to BC type map");
  auto bc_configs = handler["boundary_conditions"];
  std::map<std::string_view, std::string_view> bc_type_map;
  for (size_t i = 0; i < bc_configs.length(); i++) {
    auto bc_config = bc_configs[i];
    bc_type_map.insert(
        {bc_config.get_string("bc_name"), bc_config.get_string("bc_type")});
  }

  if (handler.has_key("external_solver_coupling_blocks")) {
    DEBUG_MSG("Create external coupling blocks");
    auto coupling_configs = handler["external_solver_coupling_blocks"];
    for (size_t i = 0; i < coupling_configs.length(); i++) {
      auto coupling_config = coupling_configs[i];
      auto coupling_type = coupling_config.get_string("type");
      auto coupling_name = coupling_config.get_string("name");
      auto coupling_loc = coupling_config.get_string("location");
      bool periodic = coupling_config.get_bool("periodic", true);
      auto coupling_values = coupling_config["values"];

      // Create coupling block
      auto t_coupling = coupling_values.get_double_array("t", {0.0});

      if (coupling_type == "FLOW") {
        auto Q_coupling = coupling_values.get_double_array("Q");

        MODEL::TimeDependentParameter q_coupling_param(t_coupling, Q_coupling,
                                                       periodic);
        if ((q_coupling_param.isconstant == false) &&
            (q_coupling_param.isperiodic == true)) {
          if ((sim_cardiac_cycle_period > 0.0) &&
              (q_coupling_param.cycle_period != sim_cardiac_cycle_period)) {
            throw std::runtime_error(
                "Inconsistent cardiac cycle period defined in "
                "TimeDependentParameter");
          } else {
            sim_cardiac_cycle_period = q_coupling_param.cycle_period;
          }
        }
        model->blocks.push_back(new MODEL::FlowReferenceBC<T>(
            q_coupling_param, static_cast<std::string>(coupling_name),
            static_cast<std::string>(coupling_loc)));
        model->block_index_map.insert(
            {static_cast<std::string>(coupling_name), block_count});
        block_count++;
        model->external_coupling_blocks.push_back(
            static_cast<std::string>(coupling_name));
      } else if (coupling_type == "PRESSURE") {
        auto P_coupling = coupling_values.get_double_array("P");
        MODEL::TimeDependentParameter p_coupling_param(t_coupling, P_coupling,
                                                       periodic);
        if ((p_coupling_param.isconstant == false) &&
            (p_coupling_param.isperiodic == true)) {
          if ((sim_cardiac_cycle_period > 0.0) &&
              (p_coupling_param.cycle_period != sim_cardiac_cycle_period)) {
            throw std::runtime_error(
                "Inconsistent cardiac cycle period defined in "
                "TimeDependentParameter");
          } else {
            sim_cardiac_cycle_period = p_coupling_param.cycle_period;
          }
        }
        model->blocks.push_back(new MODEL::PressureReferenceBC<T>(
            p_coupling_param, static_cast<std::string>(coupling_name),
            static_cast<std::string>(coupling_loc)));
        model->block_index_map.insert(
            {static_cast<std::string>(coupling_name), block_count});
        block_count++;
        model->external_coupling_blocks.push_back(
            static_cast<std::string>(coupling_name));
      } else {
        throw std::runtime_error(
            "Error. Flowsolver coupling block types should be FLOW or "
            "PRESSURE.");
      }
      DEBUG_MSG("Created coupling block " << coupling_name);

      // Determine the type of connected block
      auto connected_block = coupling_config.get_string("connected_block");
      std::string_view connected_type;
      int found_block = 0;
      if (connected_block == "ClosedLoopHeartAndPulmonary") {
        connected_type = "ClosedLoopHeartAndPulmonary";
        found_block = 1;
      } else {
        try {
          connected_type = bc_type_map.at(connected_block);
          found_block = 1;
        } catch (...) {
        }
        if (found_block == 0) {
          // Search for connected_block in the list of vessel names
          for (auto const vessel : vessel_id_map) {
            if (connected_block == vessel.second) {
              connected_type = "BloodVessel";
              found_block = 1;
              break;
            }
          }
        }
        if (found_block == 0) {
          std::cout << "Error! Could not connected type for block: "
                    << connected_block << std::endl;
          throw std::runtime_error("Terminating.");
        }
      }  // connected_block != "ClosedLoopHeartAndPulmonary"
      // Create connections
      if (coupling_loc == "inlet") {
        std::vector<std::string_view> possible_types = {
            "RESISTANCE",    "RCR",      "ClosedLoopRCR",
            "SimplifiedRCR", "CORONARY", "ClosedLoopCoronary",
            "BloodVessel"};
        if (std::find(std::begin(possible_types), std::end(possible_types),
                      connected_type) == std::end(possible_types)) {
          throw std::runtime_error(
              "Error: The specified connection type for inlet "
              "external_coupling_block is invalid.");
        }
        connections.push_back({coupling_name, connected_block});
        DEBUG_MSG("Created coupling block connection: " << coupling_name << "->"
                                                        << connected_block);
      } else if (coupling_loc == "outlet") {
        std::vector<std::string_view> possible_types = {
            "ClosedLoopRCR", "ClosedLoopHeartAndPulmonary", "BloodVessel"};
        if (std::find(std::begin(possible_types), std::end(possible_types),
                      connected_type) == std::end(possible_types)) {
          throw std::runtime_error(
              "Error: The specified connection type for outlet "
              "external_coupling_block is invalid.");
        }
        // Add connection only for closedLoopRCR and BloodVessel. Connection to
        // ClosedLoopHeartAndPulmonary will be handled in
        // ClosedLoopHeartAndPulmonary creation.
        if ((connected_type == "ClosedLoopRCR") ||
            (connected_type == "BloodVessel")) {
          connections.push_back({connected_block, coupling_name});
          DEBUG_MSG("Created coupling block connection: "
                    << connected_block << "-> " << coupling_name);
        }  // connected_type == "ClosedLoopRCR"
      }    // coupling_loc
    }      // for (size_t i = 0; i < coupling_configs.length(); i++)
  }

  std::vector<std::string_view> closed_loop_bcs;

  // Create boundary conditions
  DEBUG_MSG("Create boundary conditions");
  for (size_t i = 0; i < bc_configs.length(); i++) {
    auto bc_config = bc_configs[i];
    auto bc_type = bc_config.get_string("bc_type");
    auto bc_name = bc_config.get_string("bc_name");
    auto bc_values = bc_config["bc_values"];

    auto t = bc_values.get_double_array("t", {0.0});

    if (bc_type == "RCR") {
      T Rp = bc_values.get_double("Rp");
      T C = bc_values.get_double("C");
      T Rd = bc_values.get_double("Rd");
      T Pd = bc_values.get_double("Pd");
      model->blocks.push_back(new MODEL::WindkesselBC<T>(
          Rp = Rp, C = C, Rd = Rd, Pd = Pd, static_cast<std::string>(bc_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(bc_name), block_count});

    } else if (bc_type == "ClosedLoopRCR") {
      T Rp = bc_values.get_double("Rp");
      T C = bc_values.get_double("C");
      T Rd = bc_values.get_double("Rd");
      bool closed_loop_outlet = bc_values.get_bool("closed_loop_outlet");
      if (closed_loop_outlet == true) {
        closed_loop_bcs.push_back(bc_name);
      }
      model->blocks.push_back(new MODEL::ClosedLoopRCRBC<T>(
          Rp = Rp, C = C, Rd = Rd, closed_loop_outlet = closed_loop_outlet,
          static_cast<std::string>(bc_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(bc_name), block_count});
    } else if (bc_type == "FLOW") {
      auto Q = bc_values.get_double_array("Q");

      MODEL::TimeDependentParameter q_param(t, Q);
      if ((q_param.isconstant == false) && (q_param.isperiodic == true)) {
        if ((sim_cardiac_cycle_period > 0.0) &&
            (q_param.cycle_period != sim_cardiac_cycle_period)) {
          throw std::runtime_error(
              "Inconsistent cardiac cycle period defined in "
              "TimeDependentParameter");
        } else {
          sim_cardiac_cycle_period = q_param.cycle_period;
        }
      }
      model->blocks.push_back(new MODEL::FlowReferenceBC<T>(
          q_param, static_cast<std::string>(bc_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(bc_name), block_count});

    } else if (bc_type == "RESISTANCE") {
      auto R = bc_values.get_double_array("R");
      MODEL::TimeDependentParameter r_param(t, R);

      auto Pd = bc_values.get_double_array("Pd");
      MODEL::TimeDependentParameter pd_param(t, Pd);
      model->blocks.push_back(new MODEL::ResistanceBC<T>(
          r_param, pd_param, static_cast<std::string>(bc_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(bc_name), block_count});
    } else if (bc_type == "PRESSURE") {
      auto P = bc_values.get_double_array("P");
      MODEL::TimeDependentParameter p_param(t, P);
      if ((p_param.isconstant == false) && (p_param.isperiodic == true)) {
        if ((sim_cardiac_cycle_period > 0.0) &&
            (p_param.cycle_period != sim_cardiac_cycle_period)) {
          throw std::runtime_error(
              "Inconsistent cardiac cycle period defined in "
              "TimeDependentParameter");
        } else {
          sim_cardiac_cycle_period = p_param.cycle_period;
        }
      }
      model->blocks.push_back(new MODEL::PressureReferenceBC<T>(
          p_param, static_cast<std::string>(bc_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(bc_name), block_count});
    } else if (bc_type == "CORONARY") {
      T Ra = bc_values.get_double("Ra1");
      T Ram = bc_values.get_double("Ra2");
      T Rv = bc_values.get_double("Rv1");
      T Ca = bc_values.get_double("Ca");
      T Cim = bc_values.get_double("Cc");

      auto Pim = bc_values.get_double_array("Pim");
      MODEL::TimeDependentParameter pim_param(t, Pim);
      auto P_v = bc_values.get_double_array("P_v");
      MODEL::TimeDependentParameter pv_param(t, P_v);

      model->blocks.push_back(new MODEL::OpenLoopCoronaryBC<T>(
          Ra = Ra, Ram = Ram, Rv = Rv, Ca = Ca, Cim = Cim, pim_param, pv_param,
          static_cast<std::string>(bc_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(bc_name), block_count});
    } else if (bc_type == "ClosedLoopCoronary") {
      auto Ra = bc_values.get_double("Ra");
      auto Ram = bc_values.get_double("Ram");
      auto Rv = bc_values.get_double("Rv");
      auto Ca = bc_values.get_double("Ca");
      auto Cim = bc_values.get_double("Cim");
      auto side = bc_values.get_string("side");
      closed_loop_bcs.push_back(bc_name);
      model->blocks.push_back(new MODEL::ClosedLoopCoronaryBC<T>(
          Ra = Ra, Ram = Ram, Rv = Rv, Ca = Ca, Cim = Cim,
          static_cast<std::string>(side), static_cast<std::string>(bc_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(bc_name), block_count});
    } else {
      throw std::invalid_argument("Unknown boundary condition type");
    }
    block_count++;
    DEBUG_MSG("Created boundary condition " << bc_name);
  }

  // Create junctions
  auto junctions = handler["junctions"];
  for (size_t i = 0; i < junctions.length(); i++) {
    auto junction_config = junctions[i];
    auto j_type = junction_config.get_string("junction_type");
    auto junction_name = junction_config.get_string("junction_name");
    if ((j_type == "NORMAL_JUNCTION") || (j_type == "internal_junction")) {
      model->blocks.push_back(
          new MODEL::Junction<T>(static_cast<std::string>(junction_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(junction_name), block_count});
    } else if (j_type == "resistive_junction") {
      auto junction_values = junction_config["junction_values"];
      auto R = junction_values.get_double_array("R");
      model->blocks.push_back(new MODEL::ResistiveJunction<T>(
          R, static_cast<std::string>(junction_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(junction_name), block_count});
    } else if (j_type == "BloodVesselJunction") {
      auto junction_values = junction_config["junction_values"];
      auto R = junction_values.get_double_array("R_poiseuille");
      auto C = junction_values.get_double_array("C");
      auto L = junction_values.get_double_array("L");
      auto stenosis_coefficient =
          junction_values.get_double_array("stenosis_coefficient");
      model->blocks.push_back(new MODEL::BloodVesselJunction<T>(
          R, C, L, stenosis_coefficient,
          static_cast<std::string>(junction_name)));
      model->block_index_map.insert(
          {static_cast<std::string>(junction_name), block_count});
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
    block_count++;
    DEBUG_MSG("Created junction " << junction_name);
  }

  // Create closed-loop blocks
  bool heartpulmonary_block_present =
      false;  ///< Flag to check if heart block is present (requires different
              ///< handling)
  if (handler.has_key("closed_loop_blocks")) {
    auto closed_loop_configs = handler["closed_loop_blocks"];
    for (size_t i = 0; i < closed_loop_configs.length(); i++) {
      auto closed_loop_config = closed_loop_configs[i];
      auto closed_loop_type = closed_loop_config.get_string("closed_loop_type");
      if (closed_loop_type == "ClosedLoopHeartAndPulmonary") {
        if (heartpulmonary_block_present == false) {
          heartpulmonary_block_present = true;
          if (sim_steady_initial == true) {
            std::runtime_error(
                "ERROR: Steady initial condition is not compatible with "
                "ClosedLoopHeartAndPulmonary block.");
          }
          sim_steady_initial = false;
          std::string_view heartpulmonary_name = "CLH";
          T cycle_period =
              closed_loop_config.get_double("cardiac_cycle_period");
          if ((sim_cardiac_cycle_period > 0.0) &&
              (cycle_period != sim_cardiac_cycle_period)) {
            throw std::runtime_error(
                "Inconsistent cardiac cycle period defined in "
                "ClosedLoopHeartAndPulmonary.");
          } else {
            sim_cardiac_cycle_period = cycle_period;
          }
          auto heart_params = closed_loop_config["parameters"];
          // Convert to std::map to keep model blocks independent of simdjson
          std::map<std::string, T> param_values;
          param_values.insert(
              std::make_pair("Tsa", heart_params.get_double("Tsa")));
          param_values.insert(
              std::make_pair("tpwave", heart_params.get_double("tpwave")));
          param_values.insert(
              std::make_pair("Erv_s", heart_params.get_double("Erv_s")));
          param_values.insert(
              std::make_pair("Elv_s", heart_params.get_double("Elv_s")));
          param_values.insert(
              std::make_pair("iml", heart_params.get_double("iml")));
          param_values.insert(
              std::make_pair("imr", heart_params.get_double("imr")));
          param_values.insert(
              std::make_pair("Lra_v", heart_params.get_double("Lra_v")));
          param_values.insert(
              std::make_pair("Rra_v", heart_params.get_double("Rra_v")));
          param_values.insert(
              std::make_pair("Lrv_a", heart_params.get_double("Lrv_a")));
          param_values.insert(
              std::make_pair("Rrv_a", heart_params.get_double("Rrv_a")));
          param_values.insert(
              std::make_pair("Lla_v", heart_params.get_double("Lla_v")));
          param_values.insert(
              std::make_pair("Rla_v", heart_params.get_double("Rla_v")));
          param_values.insert(
              std::make_pair("Llv_a", heart_params.get_double("Llv_a")));
          param_values.insert(
              std::make_pair("Rlv_ao", heart_params.get_double("Rlv_ao")));
          param_values.insert(
              std::make_pair("Vrv_u", heart_params.get_double("Vrv_u")));
          param_values.insert(
              std::make_pair("Vlv_u", heart_params.get_double("Vlv_u")));
          param_values.insert(
              std::make_pair("Rpd", heart_params.get_double("Rpd")));
          param_values.insert(
              std::make_pair("Cp", heart_params.get_double("Cp")));
          param_values.insert(
              std::make_pair("Cpa", heart_params.get_double("Cpa")));
          param_values.insert(
              std::make_pair("Kxp_ra", heart_params.get_double("Kxp_ra")));
          param_values.insert(
              std::make_pair("Kxv_ra", heart_params.get_double("Kxv_ra")));
          param_values.insert(
              std::make_pair("Kxp_la", heart_params.get_double("Kxp_la")));
          param_values.insert(
              std::make_pair("Kxv_la", heart_params.get_double("Kxv_la")));
          param_values.insert(
              std::make_pair("Emax_ra", heart_params.get_double("Emax_ra")));
          param_values.insert(
              std::make_pair("Emax_la", heart_params.get_double("Emax_la")));
          param_values.insert(
              std::make_pair("Vaso_ra", heart_params.get_double("Vaso_ra")));
          param_values.insert(
              std::make_pair("Vaso_la", heart_params.get_double("Vaso_la")));
          if (param_values.size() == 27) {
            model->blocks.push_back(new MODEL::ClosedLoopHeartPulmonary<T>(
                param_values, cycle_period,
                static_cast<std::string>(heartpulmonary_name)));
            model->block_index_map.insert(
                {static_cast<std::string>(heartpulmonary_name), block_count});
            block_count++;
          } else {
            throw std::runtime_error(
                "Error. ClosedLoopHeartAndPulmonary should have 27 parameters");
          }
          // Junction at inlet to heart
          std::string_view heart_inlet_junction_name = "J_heart_inlet";
          connections.push_back(
              {heart_inlet_junction_name, heartpulmonary_name});
          model->blocks.push_back(new MODEL::Junction<T>(
              static_cast<std::string>(heart_inlet_junction_name)));
          model->block_index_map.insert(
              {static_cast<std::string>(heart_inlet_junction_name),
               block_count});
          block_count++;
          for (auto heart_inlet_elem : closed_loop_bcs) {
            connections.push_back(
                {heart_inlet_elem, heart_inlet_junction_name});
          }
          // Junction at outlet from heart
          std::string_view heart_outlet_junction_name = "J_heart_outlet";
          connections.push_back(
              {heartpulmonary_name, heart_outlet_junction_name});
          model->blocks.push_back(new MODEL::Junction<T>(
              static_cast<std::string>(heart_outlet_junction_name)));
          model->block_index_map.insert(
              {static_cast<std::string>(heart_outlet_junction_name),
               block_count});
          for (auto outlet_block :
               closed_loop_config.get_string_array("outlet_blocks")) {
            connections.push_back({heart_outlet_junction_name, outlet_block});
          }
          block_count++;
        } else {
          throw std::runtime_error(
              "Error. Only one ClosedLoopHeartAndPulmonary can be included.");
        }
      }
    }
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

  // Sanity checks
  bool model_is_valid = true;
  if (block_count != model->blocks.size()) {
    model_is_valid = false;
  }
  if (block_count != model->block_index_map.size()) {
    model_is_valid = false;
  }
  if (model_is_valid == false) {
    throw std::runtime_error(
        "Model sanity check failed: An unexpected error occured while loading "
        "the configuration. This is likely an implementation issue. Please "
        "contact a developer.");
  }

  // Setup degrees of freedom of the system
  for (auto &block : model->blocks) {
    block->setup_dofs(model->dofhandler);
  }

  // Set value of cardiac cycle period
  if (sim_cardiac_cycle_period < 0.0) {
    sim_cardiac_cycle_period =
        1.0;  // If it has not been read from config or TimeDependentParameter
              // yet, set as default value of 1.0
  }
  // Calculate time step size
  if (!sim_coupled) {
    sim_time_step_size =
        sim_cardiac_cycle_period / (T(sim_pts_per_cycle) - 1.0);
  } else {
    sim_time_step_size = sim_external_step_size / (T(sim_num_time_steps) - 1.0);
  }

  // Update block parameters that depend on DOFs and other params of other
  // blocks For example, coronary block params that depend on heart block DOFs
  for (auto &block : model->blocks) {
    block->set_model_dependent_params(*model);
  }

  // Read initial condition
  initial_state = ALGEBRA::State<T>::Zero(model->dofhandler.size());
  // Initialize blocks that have fixed initial conditions
  for (auto &block : model->blocks) {
    block->set_ICs(initial_state);
  }
  if (handler.has_key("initial_condition")) {
    auto initial_condition = handler["initial_condition"];
    // Check for pressure_all or flow_all condition.
    // This will initialize all pressure:* and flow:* variables.
    T init_p, init_q;
    bool init_p_flag = initial_condition.has_key("pressure_all");
    bool init_q_flag = initial_condition.has_key("flow_all");
    if (init_p_flag) {
      init_p = initial_condition.get_double("pressure_all");
    }
    if (init_q_flag) {
      init_q = initial_condition.get_double("flow_all");
    }

    // Loop through variables and check for initial conditions.
    for (size_t i = 0; i < model->dofhandler.size(); i++) {
      std::string var_name = model->dofhandler.variables[i];
      double default_val = 0.0;
      if ((init_p_flag == true) && ((var_name.substr(0, 9) == "pressure:") ||
                                    (var_name.substr(0, 4) == "P_c:"))) {
        default_val = init_p;
        DEBUG_MSG("pressure_all initial condition for " << var_name);
      } else if ((init_q_flag == true) && (var_name.substr(0, 5) == "flow:")) {
        default_val = init_q;
        DEBUG_MSG("flow_all initial condition for " << var_name);
      } else {
        DEBUG_MSG("No initial condition found for " << var_name);
      }
      initial_state.y[i] = initial_condition.get_double(var_name, default_val);
    }
  }
  DEBUG_MSG("[configreader] END");
}

}  // namespace IO

#endif  // SVZERODSOLVER_IO_CONFIGREADER_HPP_
