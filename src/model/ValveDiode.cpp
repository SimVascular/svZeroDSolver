// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ValveDiode.h"

#include "Model.h"

double ValveDiode::get_beta(const std::vector<double>& parameters) const {
  double rmin = parameters[global_param_ids[ParamId::RMIN]];
  double rmax = parameters[global_param_ids[ParamId::RMAX]];
  // A non-positive Rmax means an ideal closed valve: no reverse flow
  return rmax > 0.0 ? rmin / rmax : 0.0;
}

void ValveDiode::setup_model_dependent_params() {
  double rmin = model->get_parameter_value(global_param_ids[ParamId::RMIN]);
  double rmax = model->get_parameter_value(global_param_ids[ParamId::RMAX]);
  if (rmin <= 0.0 || (rmax > 0.0 && rmax <= rmin)) {
    throw std::runtime_error("ValveDiode '" + get_name() +
                             "': need Rmin > 0 and, if given, Rmax > Rmin.");
  }
}

void ValveDiode::setup_dofs(DOFHandler& dofhandler) {
  // 2 eqns, one for pressure and one for flow. The valve status feeds back
  // into nothing, so it is a post-processing variable rather than an unknown
  Block::setup_dofs_(dofhandler, 2, {}, {"valve_status"});
}

void ValveDiode::update_constant(SparseSystem& system,
                                 std::vector<double>& parameters) {
  double rmin = parameters[global_param_ids[ParamId::RMIN]];
  double beta = get_beta(parameters);

  // Eqn 0: Rmin*q_in - 0.5*(1+beta)*dp - 0.5*(1-beta)*sqrt(dp^2 + 4 eps^2)
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -0.5 * (1.0 + beta);
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 0.5 * (1.0 + beta);
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = rmin;

  // Eqn 1: q_in - q_out = 0
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
}

void ValveDiode::update_outputs(Eigen::Matrix<double, Eigen::Dynamic, 1>& y) {
  double dp = y[global_var_ids[0]] - y[global_var_ids[2]];
  double eps = model->get_parameter_value(global_param_ids[ParamId::EPS]);
  double root = sqrt(dp * dp + 4.0 * eps * eps);
  y[global_var_ids[4]] = root > 0.0 ? 0.5 * (1.0 - dp / root) : 1.0;
}

void ValveDiode::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  // Get states
  double p_in = y[global_var_ids[0]];
  double p_out = y[global_var_ids[2]];
  double dp = p_in - p_out;

  // Get parameters
  double eps = parameters[global_param_ids[ParamId::EPS]];
  double beta = get_beta(parameters);

  double root = sqrt(dp * dp + 4.0 * eps * eps);
  // Slope of the smoothed maximum, in [0, 1]. At eps = dp = 0 this picks the
  // closed branch of the subdifferential
  double s = root > 0.0 ? 0.5 * (1.0 + dp / root) : 0.0;

  // Nonlinear term
  system.C(global_eqn_ids[0]) = -0.5 * (1.0 - beta) * root;

  // Derivative of the non-linear term, d(0.5*root)/d(dp) = s - 0.5
  double droot = (1.0 - beta) * (s - 0.5);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -droot;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) = droot;
}
