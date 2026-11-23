// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file Solver.h
 * @brief Solver source file
 */

#include <string>
#include <vector>

#include "Integrator.h"
#include "Model.h"
#include "SimulationParameters.h"
#include "State.h"
#include "debug.h"

#ifndef SVZERODSOLVER_SOLVE_SOLVER_HPP_
#define SVZERODSOLVER_SOLVE_SOLVER_HPP_

/**
 * @brief Class for running 0D simulations.
 *
 * The solver solves for pressure and flow rate at the nodes of the
 * lumped-parameter system \cite pfaller22. The lumped-parameter network is
 * documented in SparseSystem. Refer to Integrator for the time discretization.
 *
 */
class Solver {
 public:
  /**
   * @brief Construct a new Solver object
   *
   * @param config Configuration handler
   */
  Solver(const nlohmann::json& config);

  /// Set up and initialize the simulation parameters and model
  void setup_initial();

  /// Set up integrator
  void setup_integrator();

  /// Run the integration
  void run_integration();

  /// Run the simulation
  void run();

  /**
   * @brief Get the full result as a csv encoded string
   *
   * @return std::string Result
   */
  std::string get_full_result() const;

  /**
   * @brief Get the result of a single DOF over time
   *
   * @param dof_name Name of the degree-of-freedom
   * @return Eigen::VectorXd Result
   */
  Eigen::VectorXd get_single_result(const std::string& dof_name) const;

  /**
   * @brief Get the result of a single DOF averaged over time
   *
   * @param dof_name Name of the degree-of-freedom
   * @return double Result
   */
  double get_single_result_avg(const std::string& dof_name) const;

  /**
   * @brief Get the time steps of the result
   *
   * @return std::vector<double> Vector of times
   */
  std::vector<double> get_times() const;

  /**
   * @brief Update the parameters of a block
   *
   * @param block_name Name of the block
   * @param new_params New parameters
   */
  void update_block_params(const std::string& block_name,
                           const std::vector<double>& new_params);

  /**
   * @brief Read the parameters of a block
   *
   * @param block_name Name of the block
   *
   * @return std::vector<double> Block parameters
   */
  std::vector<double> read_block_params(const std::string& block_name);

  /**
   * @brief Write the result to a csv file.
   *
   * @param filename
   */
  void write_result_to_csv(const std::string& filename) const;

  /**
   * @brief Get the convergence history of the whole run
   *
   * The history starts with the steps of the steady initial condition (if one
   * was computed) and continues with the steps of the transient simulation.
   * Use get_step_phases() to tell the two apart.
   *
   * @return ConvergenceHistory Per-step convergence history
   */
  ConvergenceHistory get_convergence_history() const;

  /**
   * @brief Get the number of Newton iterations used in each time step
   *
   * @return std::vector<int> Newton iterations per time step
   */
  std::vector<int> get_nonlinear_iterations() const;

  /**
   * @brief Get the largest residual at the end of each time step
   *
   * @return std::vector<double> Residual per time step
   */
  std::vector<double> get_step_residuals() const;

  /**
   * @brief Get the equation holding the largest residual of each time step
   *
   * @return std::vector<std::string> Equation name per time step
   */
  std::vector<std::string> get_step_residual_equations() const;

  /**
   * @brief Get the time at the start of each integration step
   *
   * @return std::vector<double> Time per time step
   */
  std::vector<double> get_step_times() const;

  /**
   * @brief Get the phase each step belongs to
   *
   * Either "steady_initial" for the steps of the steady initial condition or
   * "transient" for the steps of the simulation itself.
   *
   * @return std::vector<std::string> Phase per time step
   */
  std::vector<std::string> get_step_phases() const;

  /**
   * @brief Get the convergence history as a csv encoded string
   *
   * @return std::string Convergence history
   */
  std::string get_convergence_history_csv() const;

  /**
   * @brief Write the convergence history to a csv file.
   *
   * @param filename
   */
  void write_convergence_history_to_csv(const std::string& filename) const;

  /**
   * @brief Whether the configuration asks for the convergence history to be
   * written to a file
   *
   * @return bool True if `simulation_parameters.output_convergence_history` is
   * set
   */
  bool output_convergence_history() const;

 private:
  std::shared_ptr<Model> model;
  SimulationParameters simparams;
  std::vector<State> states;
  State initial_state;
  State state;
  std::vector<double> times;
  double time;
  Integrator integrator;
  ConvergenceHistory steady_history;  ///< History of the steady initialization

  void sanity_checks();

  /**
   * @brief Get indices of flow and pressure degrees-of-freedom in solution
   * vector for all vessel caps
   *
   * @return std::vector<std::pair<int, int>> Indices of flow and pressure
   * degrees-of-freedom in solution vector for all vessel caps
   */
  std::vector<std::pair<int, int>> get_vessel_caps_dof_indices();

  /**
   * @brief Check if flows and pressures for all vessel caps have converged,
   * based on cycle-to-cycle error for last two simulated cardiac cycles
   *
   * @param states_last_two_cycles Vector of solution states for last two
   * simulated cardiac cycles
   * @param vessel_caps_dof_indices Indices of flow and pressure
   * degrees-of-freedom in solution vector for all vessel caps
   *
   * @return bool True if flows and pressures for all vessel caps have converged
   */
  bool check_vessel_cap_convergence(
      const std::vector<State>& states_last_two_cycles,
      const std::vector<std::pair<int, int>>& vessel_caps_dof_indices);

  /**
   * @brief Get cycle-to-cycle errors for flow and pressure for a single vessel
   * cap
   *
   * @param states_last_two_cycles Vector of solution states for last two
   * simulated cardiac cycles
   * @param dof_indices Indices of flow and pressure degrees-of-freedom in
   * solution vector for a single vessel cap
   *
   * @return std::pair<double, double> Cycle-to-cycle errors for flow and
   * pressure
   */
  std::pair<double, double> get_cycle_to_cycle_errors_in_flow_and_pressure(
      const std::vector<State>& states_last_two_cycles,
      const std::pair<int, int>& dof_indices);
};

/**
 * @brief Name of the convergence history file belonging to a result file
 *
 * Replaces a trailing `.csv` by `_convergence.csv` so that the convergence
 * history ends up next to the result it belongs to.
 *
 * @param result_filename Name of the result file
 * @return std::string Name of the convergence history file
 */
std::string convergence_history_filename(const std::string& result_filename);

#endif
