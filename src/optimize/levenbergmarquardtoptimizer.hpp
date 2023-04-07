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

  void step(Eigen::Matrix<T, Eigen::Dynamic, 1>& alpha,
            std::vector<std::vector<T>>& y_obs,
            std::vector<std::vector<T>>& dy_obs);

 private:
  Eigen::SparseMatrix<T> jacobian;
  Eigen::Matrix<T, Eigen::Dynamic, 1> residual;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> identity;
  MODEL::Model<T>* model;
  T lambda;

  int num_obs;
  int num_params;
  int num_eqns;
  int num_vars;
  int num_dpoints;
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
  identity = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(
      num_params, num_params);
}

template <typename T>
LevenbergMarquardtOptimizer<T>::~LevenbergMarquardtOptimizer() {}

template <typename T>
void LevenbergMarquardtOptimizer<T>::step(
    Eigen::Matrix<T, Eigen::Dynamic, 1>& alpha,
    std::vector<std::vector<T>>& y_obs, std::vector<std::vector<T>>& dy_obs) {
  // Assemble gradient and residual
  for (size_t i = 0; i < num_obs; i++) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> y_i =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(num_vars);
    Eigen::Matrix<T, Eigen::Dynamic, 1> dy_i =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(num_vars);
    for (size_t k = 0; k < num_vars; k++) {
      y_i[k] = y_obs[i][k];
      dy_i[k] = dy_obs[i][k];
    }

    for (size_t j = 0; j < model->get_num_blocks(true); j++) {
      auto block = model->get_block(j);
      for (size_t l = 0; l < block->global_eqn_ids.size(); l++) {
        block->global_eqn_ids[l] += num_eqns * i;
      }
      block->update_gradient(jacobian, residual, alpha, y_i, dy_i);
      for (size_t l = 0; l < block->global_eqn_ids.size(); l++) {
        block->global_eqn_ids[l] -= num_eqns * i;
      }
    }
  }

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat =
      jacobian.transpose() * jacobian + lambda * identity;
  Eigen::Matrix<T, Eigen::Dynamic, 1> vec = jacobian.transpose() * residual;
  Eigen::Matrix<T, Eigen::Dynamic, 1> delta =
      mat.colPivHouseholderQr().solve(vec);  // or llt()

  T residual_norm = 0.0;
  for (size_t i = 0; i < num_dpoints; i++) {
    residual_norm += fabs(residual[i]);
  }
  residual_norm /= num_dpoints;
  std::cout << "residual norm: " << residual_norm << std::endl;

  T param_norm = 0.0;
  for (size_t i = 0; i < num_params; i++) {
    param_norm += fabs(delta[i]);
  }
  param_norm /= num_params;
  std::cout << "param norm: " << param_norm << std::endl;

  alpha -= delta;
}

}  // namespace OPT

#endif  // SVZERODSOLVER_OPTIMIZE_LEVENBERGMARQUARDT_HPP_
