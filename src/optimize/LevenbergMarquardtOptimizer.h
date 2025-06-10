// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file LevenbergMarquardtOptimizer.h
 * @brief opt::LevenbergMarquardtOptimizer source file
 */
#ifndef SVZERODSOLVER_OPTIMIZE_LEVENBERGMARQUARDT_HPP_
#define SVZERODSOLVER_OPTIMIZE_LEVENBERGMARQUARDT_HPP_

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Model.h"

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
 */
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
  LevenbergMarquardtOptimizer(Model* model, int num_obs, int num_params,
                              double lambda0, double tol_grad, double tol_inc,
                              int max_iter);

  /**
   * @brief Run the optimization algorithm
   *
   * @param alpha Initial parameter vector alpha
   * @param y_obs Matrix (num_obs x n) with all observations for y
   * @param dy_obs Matrix (num_obs x n) with all observations for dy
   * @return Eigen::Matrix<double, Eigen::Dynamic, 1> Optimized parameter vector
   * alpha
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> run(
      Eigen::Matrix<double, Eigen::Dynamic, 1> alpha,
      std::vector<std::vector<double>>& y_obs,
      std::vector<std::vector<double>>& dy_obs);

 private:
  Eigen::SparseMatrix<double> jacobian;
  Eigen::Matrix<double, Eigen::Dynamic, 1> residual;
  Eigen::Matrix<double, Eigen::Dynamic, 1> delta;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat;
  Eigen::Matrix<double, Eigen::Dynamic, 1> vec;
  Model* model;
  double lambda;

  int num_obs;
  int num_params;
  int num_eqns;
  int num_vars;
  int num_dpoints;

  double tol_grad;
  double tol_inc;
  int max_iter;

  void update_gradient(Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha,
                       std::vector<std::vector<double>>& y_obs,
                       std::vector<std::vector<double>>& dy_obs);

  void update_delta(bool first_step);
};

#endif  // SVZERODSOLVER_OPTIMIZE_LEVENBERGMARQUARDT_HPP_
