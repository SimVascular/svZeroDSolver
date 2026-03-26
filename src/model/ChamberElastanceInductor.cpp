// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ChamberElastanceInductor.h"

void ChamberElastanceInductor::setup_dofs(DOFHandler& dofhandler) {
  // Internal variable is chamber volume
  Block::setup_dofs_(dofhandler, 3, {"Vc"});
}

void ChamberElastanceInductor::update_constant(
    SparseSystem& system, std::vector<double>& parameters) {
  double L = parameters[global_param_ids[ParamId::IMPEDANCE]];

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

void ChamberElastanceInductor::update_time(SparseSystem& system,
                                           std::vector<double>& parameters) {
  get_elastance_values(parameters);

  // Eq 0: F[0][4] = -Elas (coefficient on Vc)
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -1 * Elas;

  // In exponential passive mode (Kxp > 0), C[0] is set in update_solution.
  // In standard linear mode, set C[0] here.
  double Kxp = parameters[global_param_ids[ParamId::KXP]];
  if (Kxp <= 0.0) {
    system.C.coeffRef(global_eqn_ids[0]) = Elas * Vrest;
  }
}

void ChamberElastanceInductor::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  // --- Exponential passive P-V (atrial mode, Kxp > 0) ---
  double Kxp = parameters[global_param_ids[ParamId::KXP]];
  if (Kxp <= 0.0) return;

  double Kxv = parameters[global_param_ids[ParamId::KXV]];
  double Vaso = parameters[global_param_ids[ParamId::VASO]];
  double Emax = parameters[global_param_ids[ParamId::EMAX]];
  double Vc = y[global_var_ids[4]];

  double exp_term = std::exp(Kxv * (Vc - Vaso));
  double psi = Kxp * (exp_term - 1.0);
  double psi_d = Kxp * Kxv * exp_term;

  system.C(global_eqn_ids[0]) = act_ * Emax * Vaso + psi * (act_ - 1.0);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      psi_d * (act_ - 1.0);
}

void ChamberElastanceInductor::get_elastance_values(
    std::vector<double>& parameters) {
  double Emax = parameters[global_param_ids[ParamId::EMAX]];
  double Emin = parameters[global_param_ids[ParamId::EMIN]];
  double Vrd = parameters[global_param_ids[ParamId::VRD]];
  double Vrs = parameters[global_param_ids[ParamId::VRS]];

  act_ = activation_func_->compute(model->time);

  double Kxp = parameters[global_param_ids[ParamId::KXP]];
  if (Kxp > 0.0) {
    Elas = act_ * Emax;
    Vrest = parameters[global_param_ids[ParamId::VASO]];
  } else {
    Vrest = (1.0 - act_) * (Vrd - Vrs) + Vrs;
    Elas = (Emax - Emin) * act_ + Emin;
  }
}

void ChamberElastanceInductor::set_activation_function(
    std::unique_ptr<ActivationFunction> af) {
  activation_func_ = std::move(af);
}
