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
 * @brief Convergence history of the non-linear iterations
 *
 * Holds one entry per time step: the time at the start of the step, the
 * number of Newton updates it took, the largest absolute residual left in
 * the solution of that step and the degree of freedom that holds it.
 */
struct ConvergenceHistory {
  std::vector<double> time;       ///< Time at the start of the step
  std::vector<int> iterations;    ///< Newton iterations of the step
  std::vector<double> residual;   ///< Largest absolute residual left
  std::vector<int> residual_dof;  ///< DOF index holding that residual

  /**
   * @brief Record one time step
   *
   * @param step_time Time at the start of the step
   * @param step_iterations Number of Newton updates performed
   * @param step_residual Largest absolute residual left
   * @param step_residual_dof DOF index holding that residual
   */
  void push_back(double step_time, int step_iterations, double step_residual,
                 int step_residual_dof) {
    time.push_back(step_time);
    iterations.push_back(step_iterations);
    residual.push_back(step_residual);
    residual_dof.push_back(step_residual_dof);
  }

  /**
   * @brief Append another history at the end of this one
   *
   * @param other History to append
   */
  void append(const ConvergenceHistory& other) {
    time.insert(time.end(), other.time.begin(), other.time.end());
    iterations.insert(iterations.end(), other.iterations.begin(),
                      other.iterations.end());
    residual.insert(residual.end(), other.residual.begin(),
                    other.residual.end());
    residual_dof.insert(residual_dof.end(), other.residual_dof.begin(),
                        other.residual_dof.end());
  }

  /**
   * @brief Discard all recorded time steps
   */
  void clear() {
    time.clear();
    iterations.clear();
    residual.clear();
    residual_dof.clear();
  }

  /**
   * @brief Number of recorded time steps
   *
   * @return Number of recorded time steps
   */
  size_t size() const { return time.size(); }
};

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
  ConvergenceHistory history;  ///< Per-step convergence history
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_af;
  Eigen::Matrix<double, Eigen::Dynamic, 1> ydot_am;
  SparseSystem system;
  Model* model{nullptr};

  /**
   * @brief Evaluate the residual of the current iterate
   *
   * Interpolates the iterates to the generalized mid-point, updates the
   * solution-dependent element contributions and assembles the residual.
   *
   * @param old_state State at the beginning of the time step
   * @param new_state Current iterate
   */
  void evaluate_residual(const State& old_state, const State& new_state);

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

  /**
   * @brief Get average number of nonlinear iterations in all step calls
   *
   * @return Average number of nonlinear iterations in all step calls
   *
   */
  double avg_nonlin_iter();

  /**
   * @brief Convergence history of all time steps performed so far
   *
   * @return Per-step convergence history
   */
  const ConvergenceHistory& get_history() const { return history; }
};

#endif  // SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
