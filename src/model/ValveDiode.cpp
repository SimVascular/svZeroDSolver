// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ValveDiode.h"

#include "Model.h"

void ValveDiode::setup_dofs(DOFHandler& dofhandler) {
  // Two equations, the diode condition and flow continuity. The valve status
  // is reported but feeds back into nothing, so it is a post-processing
  // variable: it takes no row or column in the system and is filled in
  // update_outputs. For eps = 0 it is a step function of the pressure
  // difference, which would otherwise keep a time step in which the valve
  // changes state from ever meeting the tolerance.
  Block::setup_dofs_(dofhandler, 2, {}, {"valve_status"});
}

void ValveDiode::update_constant(SparseSystem& system,
                                 std::vector<double>& parameters) {
  double r = parameters[global_param_ids[ParamId::R]];
  if (r <= 0.0) {
    throw std::runtime_error(
        "ValveDiode '" + get_name() +
        "': the open resistance R must be positive, got " + std::to_string(r) +
        ". A valve with zero open resistance does not define a flow rate "
        "when it is closed.");
  }

  // Equation 0: R * q_in - 0.5 * (p_in - p_out) - 0.5 * sqrt(...) = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -0.5;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 0.5;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = r;

  // Equation 1: q_in - q_out = 0
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

  // sqrt(dp^2 + 4 eps^2) and its derivative s = d/d(dp)
  double root = sqrt(dp * dp + 4.0 * eps * eps);
  // Slope of the smoothed maximum, in [0, 1]. For eps = 0 and dp = 0 this
  // picks the closed branch of the subdifferential.
  double s = root > 0.0 ? 0.5 * (1.0 + dp / root) : 0.0;

  // Nonlinear term
  system.C(global_eqn_ids[0]) = -0.5 * root;

  // Derivative of the non-linear term. d(0.5*root)/d(dp) = s - 0.5
  double droot = s - 0.5;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -droot;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) = droot;
}
