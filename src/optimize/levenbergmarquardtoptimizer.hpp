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
 * @file levenbergmarquardtoptimizer.hpp
 * @brief OPT::LevenbergMarquardtOptimizer source file
 */
#ifndef SVZERODSOLVER_OPTIMIZE_LEVENBERGMARQUARDT_HPP_
#define SVZERODSOLVER_OPTIMIZE_LEVENBERGMARQUARDT_HPP_

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace OPT {

/**
 * @brief Levenberg-Marquardt optimization class
 *
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class LevenbergMarquardtOptimizer {
 public:
  /**
   * @brief Construct a new LevenbergMarquardtOptimizer object
   *
   * @param model The 0D model
   * @param num_obs Number of observations in optimization
   * @param num_params Number of parameters in optimization
   * @param lambda Damping factor
   */
  LevenbergMarquardtOptimizer(MODEL::Model<T>* model, int num_obs,
                              int num_params, T lambda);
  ~LevenbergMarquardtOptimizer();

  /**
   * @brief Run the optimization algorithm
   *
   * @param alpha Initial parameter vector alpha
   * @param y_obs Matrix (num_obs x n) with all observations for y
   * @param dy_obs Matrix (num_obs x n) with all observations for dy
   * @return Eigen::Matrix<T, Eigen::Dynamic, 1> Optimized parameter vector
   * alpha
   */
  Eigen::Matrix<T, Eigen::Dynamic, 1> run(
      Eigen::Matrix<T, Eigen::Dynamic, 1> alpha,
      std::vector<std::vector<T>>& y_obs, std::vector<std::vector<T>>& dy_obs);

 private:
  Eigen::SparseMatrix<T> jacobian;
  Eigen::Matrix<T, Eigen::Dynamic, 1> residual;
  Eigen::Matrix<T, Eigen::Dynamic, 1> delta;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat;
  Eigen::Matrix<T, Eigen::Dynamic, 1> vec;
  MODEL::Model<T>* model;
  T lambda;

  int num_obs;
  int num_params;
  int num_eqns;
  int num_vars;
  int num_dpoints;

  void update_gradient(Eigen::Matrix<T, Eigen::Dynamic, 1>& alpha,
                       std::vector<std::vector<T>>& y_obs,
                       std::vector<std::vector<T>>& dy_obs);

  void update_delta(bool first_step);
};

template <typename T>
LevenbergMarquardtOptimizer<T>::LevenbergMarquardtOptimizer(
    MODEL::Model<T>* model, int num_obs, int num_params, T lambda) {
  this->model = model;
  this->num_obs = num_obs;
  this->num_params = num_params;
  this->num_eqns = model->dofhandler.get_num_equations();
  this->num_vars = model->dofhandler.get_num_variables();
  this->num_dpoints = this->num_obs * this->num_eqns;
  this->lambda = lambda;
  jacobian = Eigen::SparseMatrix<T>(num_dpoints, num_params);
  residual = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(num_dpoints);
  mat =
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(num_params, num_params);
  vec = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(num_params);
}

template <typename T>
LevenbergMarquardtOptimizer<T>::~LevenbergMarquardtOptimizer() {}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> LevenbergMarquardtOptimizer<T>::run(
    Eigen::Matrix<T, Eigen::Dynamic, 1> alpha,
    std::vector<std::vector<T>>& y_obs, std::vector<std::vector<T>>& dy_obs) {
  for (size_t nliter = 0; nliter < 100; nliter++) {
    update_gradient(alpha, y_obs, dy_obs);
    if (nliter == 0) {
      update_delta(true);
    } else {
      update_delta(false);
    }
    alpha -= delta;
    std::cout << std::setprecision(12) << "Iteration " << nliter + 1
              << " | lambda: " << lambda << " | norm inc: " << delta.norm()
              << " | norm grad: " << vec.norm() << std::endl;
  }
  return alpha;
}

template <typename T>
void LevenbergMarquardtOptimizer<T>::update_gradient(
    Eigen::Matrix<T, Eigen::Dynamic, 1>& alpha,
    std::vector<std::vector<T>>& y_obs, std::vector<std::vector<T>>& dy_obs) {
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

template <typename T>
void LevenbergMarquardtOptimizer<T>::update_delta(bool first_step) {
  // Cache old gradient vector and calulcate new one
  Eigen::Matrix<T, Eigen::Dynamic, 1> vec_old = vec;
  vec = jacobian.transpose() * residual;

  // Determine new lambda parameter from new and old gradient vector
  if (!first_step) {
    lambda *= vec.norm() / vec_old.norm();
  }

  // Determine gradient matrix
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> jacobian_sq =
      jacobian.transpose() * jacobian;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> jacobian_sq_diag =
      jacobian_sq.diagonal().asDiagonal();
  mat = jacobian_sq + lambda * jacobian_sq_diag;

  // Solve for new delta
  delta = mat.llt().solve(vec);
}

}  // namespace OPT

#endif  // SVZERODSOLVER_OPTIMIZE_LEVENBERGMARQUARDT_HPP_
