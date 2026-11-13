// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file Integrator.h
 * @brief Integrator source file
 */
#ifndef SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
#define SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_

#include <Eigen/Dense>
#include <vector>

#include "Model.h"
#include "State.h"

/**
 * @brief Generalized-alpha integrator
 *
 * This class handles the time integration scheme for solving 0D blood
 * flow system using the generalized-\f$\alpha\f$ method \cite JANSEN2000305.
 *
 * Mathematical details are available on the <a
 * href="https://simvascular.github.io/documentation/rom_simulation.html#0d-solver-theory">SimVascular
 * documentation</a>.
 */

class Integrator {
 private:
  double alpha_m{0.0};
  double alpha_f{0.0};
  double gamma{0.0};
  double time_step_size{0.0};
  double ydot_init_coeff{0.0};
  double y_coeff{0.0};
  double y_coeff_jacobian{0.0};
  double atol{0.0};
  int max_iter{0};
  bool max_iter_error_to_warning{false};
  int size{0};           ///< Length of the state vector
  int num_equations{0};  ///< Size of the non-linear system that is solved
  int n_iter{0};
  int n_nonlin_iter{0};
  std::vector<int> iter_per_step;     ///< Newton iterations per time step
  std::vector<double> res_per_step;   ///< Final max residual per time step
  std::vector<int> res_dof_per_step;  ///< DOF index of largest residual
  std::vector<double> time_per_step;  ///< Time of each step
  bool record_history{false};         ///< Record every Newton residual
  int max_line_search_steps{0};       ///< Backtracking steps (0 = off)
  double line_search_factor{1.0};     ///< Accepted growth of the correction
  std::vector<std::vector<double>> residual_history;  ///< Residual per iter
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_af;
  Eigen::Matrix<double, Eigen::Dynamic, 1> ydot_am;
  SparseSystem system;
  Model* model{nullptr};

 public:
  /**
   * @brief Construct a new Integrator object
   *
   * @param model The model to simulate
   * @param time_step_size Time step size for generalized-alpha step
   * @param rho Spectral radius for generalized-alpha step
   * @param atol Absolut tolerance for non-linear iteration termination
   * @param max_iter Maximum number of non-linear iterations
   * @param max_iter_error_to_warning If true, print a warning instead of
   * throwing error when maximum iterations is reached
   */
  Integrator(Model* model, double time_step_size, double rho, double atol,
             int max_iter, bool max_iter_error_to_warning = false);

  /**
   * @brief Construct a new Integrator object
   *
   */
  Integrator();

  /**
   * @brief Destroy the Integrator object
   *
   */
  ~Integrator();

  /**
   * @brief Delete dynamically allocated memory (in class member
   * SparseSystem<double> system).
   */
  void clean();

  /**
   * @brief Update integrator parameter and system matrices with model parameter
   * updates.
   *
   * @param time_step_size Time step size for 0D model
   */
  void update_params(double time_step_size);

  /**
   * @brief Perform a time step
   *
   * @param state Current state
   * @param time Current time
   * @return New state
   */
  State step(const State& state, double time);

 private:
  /**
   * @brief Evaluate the residual of the system at a trial state
   *
   * @param st Trial state
   * @param old_state State at the beginning of the time step
   * @param res_dof Index of the degree of freedom with the largest residual
   * @return Maximum absolute residual
   */
  double evaluate_residual(const State& st, const State& old_state,
                           int& res_dof);

 public:
  /**
   * @brief Get average number of nonlinear iterations in all step calls
   *
   * @return Average number of nonlinear iterations in all step calls
   *
   */
  double avg_nonlin_iter();

  /**
   * @brief Newton iterations used in each time step so far
   *
   * @return Number of Newton iterations of each time step
   */
  const std::vector<int>& get_iter_per_step() const { return iter_per_step; }

  /**
   * @brief Maximum residual at the end of each time step so far
   *
   * @return Largest absolute residual left at the end of each time step
   */
  const std::vector<double>& get_res_per_step() const { return res_per_step; }

  /**
   * @brief Index of the equation holding the largest residual at the end of
   * each time step so far
   *
   * @return Equation index of each time step
   */
  const std::vector<int>& get_res_dof_per_step() const {
    return res_dof_per_step;
  }

  /**
   * @brief Time at the start of each step so far
   *
   * @return Time at the start of each time step
   */
  const std::vector<double>& get_time_per_step() const { return time_per_step; }

  /**
   * @brief Record the residual after every single Newton iteration
   *
   * @param on True to record the residual of every Newton iteration
   */
  void set_record_history(bool on) { record_history = on; }

  /**
   * @brief Set the number of backtracking steps of the Newton line search
   *
   * @param n Maximum number of step halvings, 0 disables the line search
   */
  void set_max_line_search_steps(int n) { max_line_search_steps = n; }

  /**
   * @brief Set the acceptance factor of the Newton line search
   *
   * A trial step is accepted when the simplified Newton correction at the
   * trial state is at most this factor times the length of the full step.
   *
   * @param f Acceptance factor, larger values damp less
   */
  void set_line_search_factor(double f) { line_search_factor = f; }

  /**
   * @brief Residual after every Newton iteration of every time step
   *
   * @return Residual of each Newton iteration, grouped by time step
   */
  const std::vector<std::vector<double>>& get_residual_history() const {
    return residual_history;
  }

  /**
   * @brief Discard the per-step convergence history
   */
  void clear_history() {
    residual_history.clear();
    iter_per_step.clear();
    res_per_step.clear();
    res_dof_per_step.clear();
    time_per_step.clear();
  }
};

#endif  // SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
