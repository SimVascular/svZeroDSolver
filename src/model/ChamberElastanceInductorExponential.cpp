// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ChamberElastanceInductorExponential.h"

#include "Model.h"

void ChamberElastanceInductorExponential::get_elastance_values(
    std::vector<double>& parameters) {
  double Emax = parameters[global_param_ids[ParamId::EMAX]];

  act_ = activation_func_->compute(model->time);
  Elas = act_ * Emax;
}

void ChamberElastanceInductorExponential::update_time(
    SparseSystem& system, std::vector<double>& parameters) {
  get_elastance_values(parameters);

  // Eq 0: P_in - A(t)*Emax*(Vc - Vaso) - (1-A(t))*psi(Vc) = 0
  // F[0][4] = -E(t). C[0] and dC_dy[0][4] are set in update_solution.
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -Elas;
}

void ChamberElastanceInductorExponential::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  double Kxp = parameters[global_param_ids[ExponentialParamId::KXP]];
  double Kxv = parameters[global_param_ids[ExponentialParamId::KXV]];
  double Vaso = parameters[global_param_ids[ExponentialParamId::VASO]];
  double Emax = parameters[global_param_ids[ParamId::EMAX]];
  double Vc = y[global_var_ids[4]];

  double exp_term = std::exp(Kxv * (Vc - Vaso));
  double psi = Kxp * (exp_term - 1.0);
  double psi_d = Kxp * Kxv * exp_term;

  // C[0] = A(t)*Emax*Vaso + (A(t)-1)*psi(Vc)
  system.C(global_eqn_ids[0]) = act_ * Emax * Vaso + psi * (act_ - 1.0);
  // dC_dy[0][4] = (A(t)-1)*psi'(Vc)
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      psi_d * (act_ - 1.0);
}
