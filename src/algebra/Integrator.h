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
  int size{0};
  int n_iter{0};
  int n_nonlin_iter{0};
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
   */
  Integrator(Model* model, double time_step_size, double rho, double atol,
             int max_iter);

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
};

#endif  // SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
