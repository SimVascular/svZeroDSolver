// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "IdealValve.h"

void IdealValve::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {});
}

void IdealValve::update_constant(SparseSystem& system,
                                 std::vector<double>& parameters) {
  double L = parameters[global_param_ids[ParamId::IMPEDANCE]];

  // Eq 0: P_in - P_out - R*Q_in - L*dQ_in = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;
  if (L > 0.0) {
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -L;
  }

  // Eq 1: Q_in - Q_out = 0
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
}

void IdealValve::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  double p_in = y[global_var_ids[0]];
  double p_out = y[global_var_ids[2]];
  double q_in = y[global_var_ids[1]];

  double Rmin = parameters[global_param_ids[ParamId::RMIN]];
  double Rmax = parameters[global_param_ids[ParamId::RMAX]];

  // Valve closes when P_in <= P_out AND Q_in < 0 (strict on Q for init stability)
  valve_ = 1.0;
  if ((p_in <= p_out) && (q_in < 0.0)) {
    valve_ = 0.0;
  }

  double resistance = (valve_ > 0.5) ? Rmin : Rmax;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -resistance;
}
