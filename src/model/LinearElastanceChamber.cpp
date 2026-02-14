// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "LinearElastanceChamber.h"

void LinearElastanceChamber::setup_dofs(DOFHandler& dofhandler) {
  // Internal variable is chamber volume
  Block::setup_dofs_(dofhandler, 3, {"Vc"});
}

void LinearElastanceChamber::update_constant(SparseSystem& system,
                                             std::vector<double>& parameters) {
  // Eq 0: P_in - E(t)(Vc - Vrest) = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;

  // Eq 1: P_in - P_out = 0
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -1.0;

  // Eq 2: Q_in - Q_out - dVc = 0
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[3]) = -1.0;
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[4]) = -1.0;
}

void LinearElastanceChamber::update_time(SparseSystem& system,
                                         std::vector<double>& parameters) {
  get_elastance_values(parameters);

  // Eq 0: P_in - E(t)(Vc - Vrest) = P_in - E(t)*Vc + E(t)*Vrest = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -Elas;
  system.C.coeffRef(global_eqn_ids[0]) =
      Elas * parameters[global_param_ids[ParamId::VREST]];
}

void LinearElastanceChamber::get_elastance_values(
    std::vector<double>& parameters) {
  double Emax = parameters[global_param_ids[ParamId::EMAX]];
  double Epass = parameters[global_param_ids[ParamId::EPASS]];

  // Compute activation using the activation function
  double phi = activation_func_->compute(model->time);

  Elas = Epass + Emax * phi;
}

void LinearElastanceChamber::set_activation_function(
    std::unique_ptr<ActivationFunction> af) {
  activation_func_ = std::move(af);
}
