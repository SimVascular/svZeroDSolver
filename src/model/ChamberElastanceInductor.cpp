// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ChamberElastanceInductor.h"

#include "Model.h"

void ChamberElastanceInductor::setup_dofs(DOFHandler& dofhandler) {
  // Internal variable is chamber volume
  Block::setup_dofs_(dofhandler, 3, {"Vc"});
}

void ChamberElastanceInductor::update_constant(
    SparseSystem& system, std::vector<double>& parameters) {
  double L = parameters[global_param_ids[ParamId::IMPEDANCE]];

  // Eq 0: P_in - E(t)(Vc - Vrest) = 0
  // F[0][0] = 1.0 (P_in coefficient)
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  // F[0][4] = -E(t) and C[0] = E(t)*Vrest are set in update_time

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

  // Eq 0: P_in - E(t)(Vc - Vrest) = 0
  // F[0][4] = -E(t), C[0] = E(t)*Vrest
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -Elas;
  system.C.coeffRef(global_eqn_ids[0]) = Elas * Vrest;
  // F[0][0] = 1.0 (P_in coefficient) is set in update_constant
}

void ChamberElastanceInductor::get_elastance_values(
    std::vector<double>& parameters) {
  double Emax = parameters[global_param_ids[ParamId::EMAX]];
  double Emin = parameters[global_param_ids[ParamId::EMIN]];
  double Vrd = parameters[global_param_ids[ParamId::VRD]];
  double Vrs = parameters[global_param_ids[ParamId::VRS]];

  act_ = activation_func_->compute(model->time);
  Vrest = (1.0 - act_) * (Vrd - Vrs) + Vrs;
  Elas = (Emax - Emin) * act_ + Emin;
}

void ChamberElastanceInductor::set_activation_function(
    std::unique_ptr<ActivationFunction> af) {
  activation_func_ = std::move(af);
}
