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
 * The 0D residual (assuming no time-dependency in parameters) is
 *
 * \f[
 * \boldsymbol{r}(\boldsymbol{\alpha}, \boldsymbol{y}, \boldsymbol{\dot{y}}) =
 * \boldsymbol{E}(\boldsymbol{\alpha}, \boldsymbol{y}) \cdot
 * \dot{\boldsymbol{y}}+\boldsymbol{F}(\boldsymbol{\alpha}, \boldsymbol{y})
 * \cdot \boldsymbol{y}+\boldsymbol{c}(\boldsymbol{\alpha}, \boldsymbol{y}) \f]
 *
 * with solution vector \f$\boldsymbol{y} \in \mathbb{R}^{N}\f$ (flow and
 * pressure at nodes), LPN parameters \f$\boldsymbol{\alpha} \in
 * \mathbb{R}^{P}\f$, system matrices \f$\boldsymbol{E},\boldsymbol{F} \in
 * \mathbb{R}^{NxN}\f$, and system vector \f$\boldsymbol{c} \in
 * \mathbb{R}^{N}\f$.
 *
 * The least squares problem can be formulated as
 *
 * \f[
 * \min _\alpha S, \quad \mathrm { with } \quad S=\sum_i^D
 * r_i^2\left(\boldsymbol{\alpha}, y_i, \dot{y}_i\right) \f]
 *
 * with given solution vectors \f$\boldsymbol{y}\f$, \f$\boldsymbol{\dot{y}}\f$
 * at all datapoints \f$D\f$. The parameter vector is iteratively improved
 * according to
 *
 * \f[
 * \boldsymbol{\alpha}^{i+1}=\boldsymbol{\alpha}^{i}+\Delta
 * \boldsymbol{\alpha}^{i+1} \f]
 *
 * wherein the increment \f$\Delta \boldsymbol{\alpha}^{i+1} \f$ is determined
 * by solving the following system:
 *
 * \f[
 * \left[\mathbf{J}^{\mathrm{T}} \mathbf{J}+\lambda
 * \operatorname{diag}\left(\mathbf{J}^{\mathrm{T}} \mathbf{J}\right)\right]^{i}
 * \cdot \Delta \boldsymbol{\alpha}^{i+1}=-\left[\mathbf{J}^{\mathrm{T}}
 * \mathbf{r}\right]^{i}, \quad \lambda^{i}=\lambda^{i-1}
 * \cdot\left\|\left[\mathbf{J}^{\mathrm{T}} \mathbf{r}\right]^{i}\right\|_2
 * /\left\|\left[\mathbf{J}^{\mathrm{T}} \mathbf{r}\right]^{i-1}\right\|_2. \f]
 *
 * The algorithm terminates when the following tolerance thresholds are reached
 *
 * \f[
 * \left\|\left[\mathbf{J}^{\mathrm{T}}
 * \mathbf{r}\right]^{\mathrm{i}}\right\|_2<\operatorname{tol}_{\text {grad
 * }}^\alpha \text { and }\left\|\Delta
 * \boldsymbol{\alpha}^{\mathrm{i}+1}\right\|_2<\mathrm{tol}_{\text {inc
 * }}^\alpha, \f]
 *
 * The Jacobian is derived from the residual as
 *
 * \f[
 * J = \frac{\partial \boldsymbol{r}}{\partial \boldsymbol{\alpha}} =
 * \frac{\partial \mathbf{E}}{\partial \boldsymbol{\alpha}} \cdot
 * \dot{\mathbf{y}}+\frac{\partial \mathbf{F}}{\partial \boldsymbol{\alpha}}
 * \cdot \mathbf{y}+\frac{\partial \mathbf{c}}{\partial \boldsymbol{\alpha}} \f]
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
   * @param lambda0 Initial damping factor
   * @param tol_grad Gradient tolerance
   * @param tol_inc Parameter increment tolerance
   * @param max_iter Maximum iterations
   */
  LevenbergMarquardtOptimizer(MODEL::Model<T>* model, int num_obs,
                              int num_params, T lambda0, T tol_grad, T tol_inc,
                              int max_iter);

  /**
   * @brief Destroy the LevenbergMarquardtOptimizer object
   *
   */
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

  T tol_grad;
  T tol_inc;
  int max_iter;

  void update_gradient(Eigen::Matrix<T, Eigen::Dynamic, 1>& alpha,
                       std::vector<std::vector<T>>& y_obs,
                       std::vector<std::vector<T>>& dy_obs);

  void update_delta(bool first_step);
};

template <typename T>
LevenbergMarquardtOptimizer<T>::LevenbergMarquardtOptimizer(
    MODEL::Model<T>* model, int num_obs, int num_params, T lambda0, T tol_grad,
    T tol_inc, int max_iter) {
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
  for (size_t i = 0; i < max_iter; i++) {
    update_gradient(alpha, y_obs, dy_obs);
    if (i == 0) {
      update_delta(true);
    } else {
      update_delta(false);
    }
    alpha -= delta;
    T norm_grad = vec.norm();
    T norm_inc = delta.norm();
    std::cout << std::setprecision(12) << "Iteration " << i + 1
              << " | lambda: " << lambda << " | norm inc: " << norm_inc
              << " | norm grad: " << norm_grad << std::endl;
    if ((norm_grad < tol_grad) && (norm_inc < tol_inc)) {
      break;
    }
    if (i >= max_iter - 1) {
      throw std::runtime_error("Maximum number of iterations reached");
    }
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
