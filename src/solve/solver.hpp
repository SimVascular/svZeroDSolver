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
 * @file solver.hpp
 * @brief SOLVE::Solver source file
 */

#include "algebra/integrator.hpp"
#include "algebra/state.hpp"
#include "helpers/debug.hpp"
#include "helpers/endswith.hpp"
#include "io/configreader.hpp"
#include "io/csvwriter.hpp"
#include "model/model.hpp"

#ifndef SVZERODSOLVER_SOLVE_SOLVER_HPP_
#define SVZERODSOLVER_SOLVE_SOLVER_HPP_

namespace SOLVE {

/**
 * @brief Class for running 0D simulations.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class Solver {
 public:
  /**
   * @brief Construct a new Solver object
   *
   * @param handler Configuration handler
   */
  Solver(const nlohmann::json& config);

  /**
   * @brief Destroy the Solver object
   *
   */
  ~Solver();

  /**
   * @brief Run the simulation
   *
   */
  void run();

  /**
   * @brief Get the full result as a csv encoded string
   *
   * @return std::string Result
   */
  std::string get_full_result();

  /**
   * @brief Get the result of a single DOF over time
   *
   * @param dof_name Name of the degree-of-freedom
   * @return Eigen::VectorXd Result
   */
  Eigen::VectorXd get_single_result(std::string dof_name);

  /**
   * @brief Get the result of a single DOF averaged over time
   *
   * @param dof_name Name of the degree-of-freedom
   * @return T Result
   */
  T get_single_result_avg(std::string dof_name);

  /**
   * @brief Get the time steps of the result
   *
   * @return std::vector<T>
   */
  std::vector<T> get_times();

  /**
   * @brief Update the parameters of a block
   *
   * @param block_name Name of the block
   * @param new_params New parameters
   */
  void update_block_params(std::string block_name, std::vector<T> new_params);

  /**
   * @brief Write the result to a csv file.
   *
   * @param filename
   */
  void write_result_to_csv(std::string filename);

 private:
  MODEL::Model<T> model;
  IO::SimulationParameters<T> simparams;
  std::vector<ALGEBRA::State<T>> states;
  std::vector<T> times;
  ALGEBRA::State<T> inital_state;

  void sanity_checks();
};

template <typename T>
Solver<T>::Solver(const nlohmann::json& config) {
  DEBUG_MSG("Read simulation parameters");
  simparams = IO::load_simulation_params<T>(config);
  DEBUG_MSG("Load model");
  model = MODEL::Model<T>();
  IO::load_simulation_model<T>(config, model);
  DEBUG_MSG("Load initial condition");
  inital_state = IO::load_initial_condition<T>(config, model);

  DEBUG_MSG("Cardiac cycle period " << model.cardiac_cycle_period);

  // Calculate time step size
  if (!simparams.sim_coupled) {
    simparams.sim_time_step_size =
        model.cardiac_cycle_period / (T(simparams.sim_pts_per_cycle) - 1.0);
  } else {
    simparams.sim_time_step_size = simparams.sim_external_step_size /
                                   (T(simparams.sim_num_time_steps) - 1.0);
  }
  sanity_checks();
}

template <typename T>
Solver<T>::~Solver() {}

template <typename T>
void Solver<T>::run() {
  auto state = inital_state;

  // Create steady initial
  if (simparams.sim_steady_initial) {
    DEBUG_MSG("Calculate steady initial condition");
    T time_step_size_steady = model.cardiac_cycle_period / 10.0;
    model.to_steady();
    ALGEBRA::Integrator<T> integrator_steady(&model, time_step_size_steady, 0.1,
                                             simparams.sim_abs_tol,
                                             simparams.sim_nliter);
    for (int i = 0; i < 31; i++) {
      state = integrator_steady.step(state, time_step_size_steady * T(i));
    }
    model.to_unsteady();
  }

  // Set-up integrator
  DEBUG_MSG("Setup time integration");
  ALGEBRA::Integrator<T> integrator(&model, simparams.sim_time_step_size, 0.1,
                                    simparams.sim_abs_tol,
                                    simparams.sim_nliter);

  // Initialize loop
  states = std::vector<ALGEBRA::State<T>>();
  times = std::vector<T>();
  if (simparams.output_all_cycles) {
    int num_states =
        simparams.sim_num_time_steps / simparams.output_interval + 1;
    states.reserve(num_states);
    times.reserve(num_states);
  } else {
    int num_states =
        simparams.sim_pts_per_cycle / simparams.output_interval + 1;
    states.reserve(num_states);
    times.reserve(num_states);
  }
  T time = 0.0;

  // Run integrator
  DEBUG_MSG("Run time integration");
  int interval_counter = 0;
  int start_last_cycle =
      simparams.sim_num_time_steps - simparams.sim_pts_per_cycle;
  if (simparams.output_all_cycles || (0 >= start_last_cycle)) {
    times.push_back(time);
    states.push_back(std::move(state));
  }
  for (int i = 1; i < simparams.sim_num_time_steps; i++) {
    state = integrator.step(state, time);
    interval_counter += 1;
    time = simparams.sim_time_step_size * T(i);
    if ((interval_counter == simparams.output_interval) ||
        (!simparams.output_all_cycles && (i == start_last_cycle))) {
      if (simparams.output_all_cycles || (i >= start_last_cycle)) {
        times.push_back(time);
        states.push_back(std::move(state));
      }
      interval_counter = 0;
    }
  }

  // Make times start from 0
  if (!simparams.output_all_cycles) {
    T start_time = times[0];
    for (auto& time : times) {
      time -= start_time;
    }
  }
}

template <typename T>
std::vector<T> Solver<T>::get_times() {
  return times;
}

template <typename T>
std::string Solver<T>::get_full_result() {
  std::string output;
  if (simparams.output_variable_based) {
    output =
        IO::to_variable_csv<T>(times, states, model, simparams.output_mean_only,
                               simparams.output_derivative);
  } else {
    output =
        IO::to_vessel_csv<T>(times, states, model, simparams.output_mean_only,
                             simparams.output_derivative);
  }
  return output;
}

template <typename T>
Eigen::VectorXd Solver<T>::get_single_result(std::string dof_name) {
  int dof_index = model.dofhandler.get_variable_index(dof_name);
  int num_states = states.size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(num_states);
  for (size_t i = 0; i < num_states; i++) result[i] = states[i].y[dof_index];
  return result;
}

template <typename T>
T Solver<T>::get_single_result_avg(std::string dof_name) {
  int dof_index = model.dofhandler.get_variable_index(dof_name);
  int num_states = states.size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(num_states);
  for (size_t i = 0; i < num_states; i++) result[i] = states[i].y[dof_index];
  return result.mean();
}

template <typename T>
void Solver<T>::update_block_params(std::string block_name,
                                    std::vector<T> new_params) {
  auto block = model.get_block(block_name);
  if (new_params.size() != block->global_param_ids.size()) {
    throw std::runtime_error(
        "Parameter update failed! Number of provided parameters does not match "
        "with block parameters.");
  }
  for (size_t i = 0; i < new_params.size(); i++) {
    model.get_parameter(block->global_param_ids[i])->update(new_params[i]);
  }
}

template <typename T>
void Solver<T>::sanity_checks() {
  // Check that steady initial is not used with ClosedLoopHeartAndPulmonary
  if ((simparams.sim_steady_initial == true) &&
      (model.get_block("CLH") != nullptr)) {
    std::runtime_error(
        "ERROR: Steady initial condition is not compatible with "
        "ClosedLoopHeartAndPulmonary block.");
  }
}

template <typename T>
void Solver<T>::write_result_to_csv(std::string filename) {
  DEBUG_MSG("Write output");
  std::ofstream ofs(filename);
  ofs << get_full_result();
  ofs.close();
}

}  // namespace SOLVE

#endif  // SVZERODSOLVER_SOLVE_SOLVER_HPP_
