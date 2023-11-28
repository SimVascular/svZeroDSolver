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
 * @file SimulationParameters.h
 * @brief Source file to read simulation configuration
 */
#ifndef SVZERODSOLVER_SIMULATIONPARAMETERS_HPP_
#define SVZERODSOLVER_SIMULATIONPARAMETERS_HPP_

#include <list>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

#include "Model.h"
#include "State.h"
#include "debug.h"

/**
 * @brief Simulation parameters
 *
 */
struct SimulationParameters {
  // Negative value indicates this has not
  // been read from config file yet.
  double sim_time_step_size{0.0};  ///< Simulation time step size
  double sim_abs_tol{0.0};         ///< Absolute tolerance for simulation

  int sim_num_cycles{0};      ///< Number of cardiac cycles to simulate
  int sim_pts_per_cycle{0};   ///< Number of time steps per cardiac cycle
  int sim_num_time_steps{0};  ///< Total number of time steps
  int sim_nliter{0};  ///< Maximum number of non-linear iterations in time
                      ///< integration
  double sim_rho_infty{0.0};  ///< Spectral radius of generalized-alpha
  int output_interval{0};     ///< Interval of writing output

  bool sim_steady_initial{0};  ///< Start from steady solution
  bool output_variable_based{
      false};  ///< Output variable based instead of vessel based
  bool output_mean_only{false};   ///< Output only the mean value
  bool output_derivative{false};  ///< Output derivatives
  bool output_all_cycles{false};  ///< Output all cardiac cycles

  bool sim_coupled{
      false};  ///< Running 0D simulation coupled with external solver
  double sim_external_step_size{0.0};  ///< Step size of external solver if
                                       ///< running coupled
};

/**
 * @brief Generate a new block
 *
 * @param model The model that teh block is added to
 * @param config The json configuration
 * @param block_type The type of block
 * @param name The name of the block
 * @param internal Is this an internal block? This is relevane for the calibrator
 * @param periodic Is this block periodic with the cardiac cycle? This is relevant for coupling with external solvers
 * @return int The block count
 */
int generate_block(Model& model, const nlohmann::json& config,
                   const std::string& block_type, const std::string_view& name,
                   bool internal = false, bool periodic = true);

/**
 * @brief Load initial conditions from a configuration
 *
 * @param config The json configuration
 * @param model The model
 * @return State Initial configuration for the model
 */
State load_initial_condition(const nlohmann::json& config, Model& model);

/**
 * @brief Load the simulation parameters from a json configuration
 *
 * @param config The json configuration
 * @return SimulationParameters Simulation parameters read from configuration
 */
SimulationParameters load_simulation_params(const nlohmann::json& config);

/**
 * @brief Load model from a configuration
 *
 * @param config The json configuration
 * @param model The 0D model
 * @
 */
void load_simulation_model(const nlohmann::json& config, Model& model);

void validate_input(const nlohmann::json& config);

void create_vessels(
    Model& model,
    std::vector<std::tuple<std::string, std::string>>& connections,
    const nlohmann::json& config, std::map<int, std::string>& vessel_id_map);

void create_external_coupling(
    Model& model,
    std::vector<std::tuple<std::string, std::string>>& connections,
    const nlohmann::json& config, std::map<int, std::string>& vessel_id_map,
    std::map<std::string, std::string>& bc_type_map);

void create_boundary_conditions(Model& model, const nlohmann::json& config,
                                std::map<std::string, std::string>& bc_type_map,
                                std::vector<std::string>& closed_loop_bcs);

void create_junctions(
    Model& model,
    std::vector<std::tuple<std::string, std::string>>& connections,
    const nlohmann::json& config, std::map<int, std::string>& vessel_id_map);

void create_closed_loop(
    Model& model,
    std::vector<std::tuple<std::string, std::string>>& connections,
    const nlohmann::json& config, std::vector<std::string>& closed_loop_bcs);

#endif
