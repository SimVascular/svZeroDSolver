// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "SmoothValveL.h"

void SmoothValveL::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {});
}

void SmoothValveL::update_constant(SparseSystem& system,
                                   std::vector<double>& parameters) {
  double L = parameters[global_param_ids[ParamId::INDUCTANCE]];

  // Eq 0: P_in - P_out - R_eff*Q_in - L*dQ_in = 0
  // Constant: P_in, P_out, L*dQ_in. R_eff*Q_in goes in C.
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -L;

  // Eq 1: Q_in - Q_out = 0
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
}

void SmoothValveL::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  double p_in = y[global_var_ids[0]];
  double p_out = y[global_var_ids[2]];
  double q_in = y[global_var_ids[1]];

  double Rmin = parameters[global_param_ids[ParamId::RMIN]];
  double Rmax = parameters[global_param_ids[ParamId::RMAX]];
  double k = parameters[global_param_ids[ParamId::STEEPNESS]];

  // s = 0.5*(1+tanh(k*ΔP)): 1 when open, 0 when closed
  double s = 0.5 * (1.0 + tanh(k * (p_in - p_out)));

  // R_eff = Rmin*s + Rmax*(1-s)
  double R_eff = Rmin * s + Rmax * (1.0 - s);

  // C[0] = -R_eff * Q_in (nonlinear product)
  system.C(global_eqn_ids[0]) = -R_eff * q_in;

  // dC/dQ_in = -R_eff (treat sigmoid as constant during linearization)
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -R_eff;
}
