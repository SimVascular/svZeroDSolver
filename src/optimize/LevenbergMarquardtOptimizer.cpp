// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "LevenbergMarquardtOptimizer.h"

#include <iomanip>

LevenbergMarquardtOptimizer::LevenbergMarquardtOptimizer(
    Model* model, int num_obs, int num_params, double lambda0, double tol_grad,
    double tol_inc, int max_iter) {
  this->model = model;
  this->num_obs = num_obs;
  this->num_params = num_params;
  this->num_eqns = model->dofhandler.get_num_equations();
  this->num_vars = model->dofhandler.get_num_variables();
  this->num_dpoints = this->num_obs * this->num_eqns;
  this->lambda = lambda0;
  this->tol_grad = tol_grad;
  this->tol_inc = tol_inc;
  this->max_iter = max_iter;

  jacobian = Eigen::SparseMatrix<double>(num_dpoints, num_params);
  residual = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(num_dpoints);
  mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(num_params,
                                                              num_params);
  vec = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(num_params);

  // Set up the solver system so that its residual assembly can be reused for a
  // single observation at a time.
  system = SparseSystem(num_vars);
  system.reserve(model);
}

Eigen::Matrix<double, Eigen::Dynamic, 1> LevenbergMarquardtOptimizer::run(
    Eigen::Matrix<double, Eigen::Dynamic, 1> alpha,
    std::vector<std::vector<double>>& y_obs,
    std::vector<std::vector<double>>& dy_obs) {
  for (size_t i = 0; i < max_iter; i++) {
    update_gradient(alpha, y_obs, dy_obs);

    if (i == 0) {
      update_delta(true);
    } else {
      update_delta(false);
    }

    alpha -= delta;
    double norm_grad = vec.norm();
    double norm_inc = delta.norm();
    std::cout << std::setprecision(1) << std::scientific << "Iteration "
              << i + 1 << " | lambda: " << lambda << " | norm inc: " << norm_inc
              << " | norm grad: " << norm_grad << std::endl;
    if ((norm_grad < tol_grad) && (norm_inc < tol_inc)) {
      break;
    }
    if (i >= max_iter - 1) {
      std::cout << "Maximum number of iterations reached" << std::endl;
      break;
    }
  }
  return alpha;
}

void LevenbergMarquardtOptimizer::update_gradient(
    Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha,
    std::vector<std::vector<double>>& y_obs,
    std::vector<std::vector<double>>& dy_obs) {
  // Set jacobian and residual to zero
  jacobian.setZero();
  residual.setZero();

  // Push the current parameter vector into the model so that the solver's
  // residual assembly uses it.
  for (size_t k = 0; k < num_params; k++) {
    model->update_parameter_value(k, alpha[k]);
  }

  // Assemble the parameter-dependent constant system contributions (E and F)
  // once for the current parameters.
  model->update_constant(system);

  Eigen::Matrix<double, Eigen::Dynamic, 1> y_dpoint(num_vars);
  Eigen::Matrix<double, Eigen::Dynamic, 1> dy_dpoint(num_vars);

  // Assemble gradient and residual
  for (size_t i = 0; i < num_obs; i++) {
    // Copy the observation for this datapoint into Eigen vectors.
    for (size_t k = 0; k < num_vars; k++) {
      y_dpoint[k] = y_obs[i][k];
      dy_dpoint[k] = dy_obs[i][k];
    }

    // Reuse the residual of the solver: r = -C - E*ydot - F*y. This is the
    // same residual the solver assembles, so it does not need to be redefined
    // here.
    model->update_solution(system, y_dpoint, dy_dpoint);
    system.update_residual(y_dpoint, dy_dpoint);
    residual.segment(num_eqns * i, num_eqns) = system.residual;

    // Assemble the Jacobian of the residual with respect to the parameters.
    for (size_t j = 0; j < model->get_num_blocks(true); j++) {
      auto block = model->get_block(j);
      // Blocks without calibration parameters do not contribute to the
      // Jacobian (e.g. a normal junction).
      if (block->global_param_ids.empty()) {
        continue;
      }
      for (size_t l = 0; l < block->global_eqn_ids.size(); l++) {
        block->global_eqn_ids[l] += num_eqns * i;
      }
      block->update_gradient(jacobian, alpha, y_obs[i], dy_obs[i]);
      for (size_t l = 0; l < block->global_eqn_ids.size(); l++) {
        block->global_eqn_ids[l] -= num_eqns * i;
      }
    }
  }
}

void LevenbergMarquardtOptimizer::update_delta(bool first_step) {
  // Cache old gradient vector and calulcate new one
  Eigen::Matrix<double, Eigen::Dynamic, 1> vec_old = vec;
  vec = jacobian.transpose() * residual;

  // Determine new lambda parameter from new and old gradient vector
  if (!first_step) {
    lambda *= vec.norm() / vec_old.norm();
  }

  // Determine gradient matrix
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> jacobian_sq =
      jacobian.transpose() * jacobian;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> jacobian_sq_diag =
      jacobian_sq.diagonal().asDiagonal();
  mat = jacobian_sq + lambda * jacobian_sq_diag;

  // Solve for new delta
  delta = mat.llt().solve(vec);
}
