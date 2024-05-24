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
 * @file Solver.h
 * @brief Solver source file
 */

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

 private:
  std::shared_ptr<Model> model;
  SimulationParameters simparams;
  std::vector<State> states;
  std::vector<double> times;
  State initial_state;

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

#endif
