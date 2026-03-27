// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ChamberElastanceInductorExponential.h"

#include "Model.h"

void ChamberElastanceInductorExponential::get_elastance_values(
    std::vector<double>& parameters) {
  double Emax = parameters[global_param_ids[ExponentialParamId::EXP_EMAX]];

  act_ = activation_func_->compute(model->time);
  Elas = act_ * Emax;
}

void ChamberElastanceInductorExponential::update_constant(
    SparseSystem& system, std::vector<double>& parameters) {
  double L = parameters[global_param_ids[ExponentialParamId::EXP_IMPEDANCE]];

  // Eq 0: P_in - E(t)(Vc - Vrest) = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;

  // Eq 1: P_in - P_out - L*dQ_out = 0
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -1.0;
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -L;

  // Eq 2: Q_in - Q_out - dVc = 0
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[3]) = -1.0;
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[4]) = -1.0;
}

void ChamberElastanceInductorExponential::update_time(
    SparseSystem& system, std::vector<double>& parameters) {
  get_elastance_values(parameters);

  // Eq 0: F[0][4] = -E(t). C[0] is set in update_solution (nonlinear).
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -Elas;
}

void ChamberElastanceInductorExponential::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  double Kxp = parameters[global_param_ids[ExponentialParamId::KXP]];
  double Kxv = parameters[global_param_ids[ExponentialParamId::KXV]];
  double Vaso = parameters[global_param_ids[ExponentialParamId::VASO]];
  double Emax = parameters[global_param_ids[ExponentialParamId::EXP_EMAX]];
  double Vc = y[global_var_ids[4]];

  double exp_term = std::exp(Kxv * (Vc - Vaso));
  double psi = Kxp * (exp_term - 1.0);
  double psi_d = Kxp * Kxv * exp_term;

  system.C(global_eqn_ids[0]) = act_ * Emax * Vaso + psi * (act_ - 1.0);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      psi_d * (act_ - 1.0);
}
