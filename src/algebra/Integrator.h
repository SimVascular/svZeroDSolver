// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file Integrator.h
 * @brief Integrator source file
 */
#ifndef SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
#define SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_

#include <Eigen/Dense>

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
 protected:
  double alpha_m{0.0};                    ///< Coefficient alpha_m
  double alpha_f{0.0};                    ///< Coefficient alpha_f
  double gamma{0.0};                      ///< Coefficient gamma
  double time_step_size{0.0};             ///< Time step size
  double ydot_init_coeff{0.0};            ///< Predictor coefficient of ydot
  double y_coeff{0.0};                    ///< gamma * time_step_size
  double y_coeff_jacobian{0.0};           ///< alpha_f * y_coeff
  double atol{0.0};                       ///< Absolute residual tolerance
  int max_iter{0};                        ///< Maximum nonlinear iterations
  bool max_iter_error_to_warning{false};  ///< Warn instead of throwing
  int size{0};                            ///< Size of the system
  int n_iter{0};                          ///< Number of step calls
  int n_nonlin_iter{0};                   ///< Number of nonlinear iterations
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_af;     ///< y at alpha_f
  Eigen::Matrix<double, Eigen::Dynamic, 1> ydot_am;  ///< ydot at alpha_m
  SparseSystem system;                               ///< System of equations
  Model* model{nullptr};                             ///< Model to integrate

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
  virtual ~Integrator();

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
  virtual State step(const State& state, double time);

  /**
   * @brief Get average number of nonlinear iterations in all step calls
   *
   * @return Average number of nonlinear iterations in all step calls
   *
   */
  double avg_nonlin_iter();
};

/**
 * @brief Integrator for stiff active mechanics
 *
 * The generalized-\f$\alpha\f$ update with spectral radius 0 (L-stable),
 * augmented with a backtracking line search on the residual norm (damped
 * Newton) and adaptive sub-stepping: a step that does not converge is retried
 * with half the step size, recursively. The force-velocity term of the active
 * stress law of \cite genet23 needs it; plain generalized-\f$\alpha\f$ does
 * not converge. Selected with the simulation parameter `"integrator":
 * "stiff"`. Blocks with persistent memory (e.g. ImpedanceBC) see one accepted
 * step of the nominal size.
 */
class ConsistentStiffIntegrator : public Integrator {
 public:
  /**
   * @brief Construct a new stiff integrator
   *
   * @param model The model to simulate
   * @param time_step_size Nominal time step size
   * @param atol Absolute tolerance for non-linear iteration termination
   * @param max_iter Maximum number of non-linear iterations per sub-step
   */
  ConsistentStiffIntegrator(Model* model, double time_step_size, double atol,
                            int max_iter);

  /**
   * @brief Perform a time step, sub-stepping if needed
   *
   * @param state Current state
   * @param time Current time
   * @return New state
   */
  State step(const State& state, double time) override;

 private:
  int max_subdiv{9};      ///< Maximum number of step halvings
  int max_backtrack{12};  ///< Maximum line-search iterations

  /**
   * @brief Residual norm at the generalized mid-point of a candidate state
   *
   * @param old_state State at the beginning of the step
   * @param new_state Candidate state at the end of the step
   * @return Infinity norm of the residual
   */
  double residual_norm(const State& old_state, const State& new_state);

  /**
   * @brief One damped-Newton step of size `dt` starting at `time`
   *
   * @param old_state State at the beginning of the sub-step
   * @param new_state State at the end of the sub-step (output)
   * @param dt Sub-step size
   * @param time Time at the beginning of the sub-step
   * @return true on convergence
   */
  bool damped_step(const State& old_state, State& new_state, double dt,
                   double time);
};

#endif  // SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
