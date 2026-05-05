// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "LevenbergMarquardtOptimizer.h"

#include <iomanip>

LevenbergMarquardtOptimizer::LevenbergMarquardtOptimizer(
    Model* model, int num_obs, int num_params,
    const std::vector<int>& active_param_ids, double lambda0, double tol_grad,
    double tol_inc, int max_iter)
    // Size the assembly system with max(num_eqns, num_vars). The calibrator
    // has no boundary-condition blocks so it can have num_vars > num_eqns;
    // sizing with the larger ensures coeffRef writes (eqn_id < num_eqns,
    // var_id < num_vars) stay in bounds.
    : system(std::max(model->dofhandler.get_num_equations(),
                      model->dofhandler.get_num_variables())) {
  this->model = model;
  this->num_obs = num_obs;
  this->num_params = num_params;
  this->active_param_ids = active_param_ids;
  this->num_active = static_cast<int>(active_param_ids.size());
  this->num_eqns = model->dofhandler.get_num_equations();
  this->num_vars = model->dofhandler.get_num_variables();
  this->num_dpoints = this->num_obs * this->num_eqns;
  this->lambda = lambda0;
  this->tol_grad = tol_grad;
  this->tol_inc = tol_inc;
  this->max_iter = max_iter;

  jacobian = Eigen::SparseMatrix<double>(num_dpoints, num_params);
  residual = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(num_dpoints);
  mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(num_active,
                                                              num_active);
  vec = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(num_active);

  // Establish the sparsity pattern of E, F, dC/dy, dC/dydot using dummy
  // values; subsequent calls to update_constant / update_solution overwrite
  // the entries via coeffRef without changing the pattern.
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

    for (int k = 0; k < num_active; k++) {
      alpha[active_param_ids[k]] -= delta[k];
    }
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
  jacobian.setZero();
  residual.setZero();

  // Sync the LM parameter vector into the model's parameter storage so the
  // solver assembly (update_constant / update_solution below) sees the same
  // values that the per-block Jacobians do.
  for (int p = 0; p < num_params; p++) {
    model->update_parameter_value(p, alpha[p]);
  }

  // E and F are constant in the parameters within one LM iteration; refresh
  // them once before walking the observations. ``update_time`` is skipped
  // because it would overwrite ``parameter_values`` from the model's
  // ``Parameter`` objects, undoing the alpha sync above.
  model->update_constant(system);

  // The system is sized to max(num_eqns, num_vars); pad y and dy so that
  // matrix-vector products E*dy and F*y are well-formed.
  int n = static_cast<int>(system.C.size());
  Eigen::Matrix<double, Eigen::Dynamic, 1> y =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  Eigen::Matrix<double, Eigen::Dynamic, 1> dy =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  for (size_t i = 0; i < num_obs; i++) {
    // Copy the observation into Eigen vectors so the forward-solver assembly
    // (update_solution) can reuse its existing API.
    for (int k = 0; k < num_vars; k++) {
      y[k] = y_obs[i][k];
      dy[k] = dy_obs[i][k];
    }
    model->update_solution(system, y, dy);

    // r_i = c + E * dy + F * y, identical to the residual the forward solver
    // computes (up to a sign that cancels in the LM normal equations). Only
    // the first num_eqns rows of the assembly are populated by vessel and
    // junction blocks; ignore the trailing padding that exists when the
    // calibrator's model has more variables than equations.
    Eigen::Matrix<double, Eigen::Dynamic, 1> r = system.C;
    r.noalias() += system.E * dy;
    r.noalias() += system.F * y;
    residual.segment(i * num_eqns, num_eqns) = r.head(num_eqns);

    // Each block contributes its own columns of the parameter Jacobian.
    for (size_t j = 0; j < model->get_num_blocks(true); j++) {
      auto block = model->get_block(j);
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
  // Build a reduced Jacobian containing only columns of active parameters.
  // Inactive parameters are held constant, so they do not appear in the
  // reduced normal equations.
  Eigen::SparseMatrix<double> jac_active(num_dpoints, num_active);
  std::vector<Eigen::Triplet<double>> triplets;
  for (int k = 0; k < num_active; k++) {
    int col = active_param_ids[k];
    for (Eigen::SparseMatrix<double>::InnerIterator it(jacobian, col); it;
         ++it) {
      triplets.emplace_back(it.row(), k, it.value());
    }
  }
  jac_active.setFromTriplets(triplets.begin(), triplets.end());

  // Cache old gradient vector and calulcate new one
  Eigen::Matrix<double, Eigen::Dynamic, 1> vec_old = vec;
  vec = jac_active.transpose() * residual;

  // Determine new lambda parameter from new and old gradient vector
  if (!first_step) {
    lambda *= vec.norm() / vec_old.norm();
  }

  // Determine gradient matrix
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> jacobian_sq =
      jac_active.transpose() * jac_active;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> jacobian_sq_diag =
      jacobian_sq.diagonal().asDiagonal();
  mat = jacobian_sq + lambda * jacobian_sq_diag;

  // Solve for new delta
  delta = mat.llt().solve(vec);
}
