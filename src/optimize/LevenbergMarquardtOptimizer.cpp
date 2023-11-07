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

  // Assemble gradient and residual
  for (size_t i = 0; i < num_obs; i++) {
    for (size_t j = 0; j < model->get_num_blocks(true); j++) {
      auto block = model->get_block(j);
      for (size_t l = 0; l < block->global_eqn_ids.size(); l++) {
        block->global_eqn_ids[l] += num_eqns * i;
      }
      block->update_gradient(jacobian, residual, alpha, y_obs[i], dy_obs[i]);
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
