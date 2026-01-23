// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "PiecewiseCosineChamber.h"

void PiecewiseCosineChamber::setup_dofs(DOFHandler& dofhandler) {
  // Internal variable is chamber volume
  Block::setup_dofs_(dofhandler, 3, {"Vc"});
}

void PiecewiseCosineChamber::update_constant(SparseSystem& system,
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

void PiecewiseCosineChamber::update_time(SparseSystem& system,
                                         std::vector<double>& parameters) {
  get_elastance_values(parameters);

  // Eq 0: P_in - E(t)(Vc - Vrest) = P_in - E(t)*Vc + E(t)*Vrest = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -1 * Elas;
  system.C.coeffRef(global_eqn_ids[0]) =
      Elas * parameters[global_param_ids[ParamId::VREST]];
}

void PiecewiseCosineChamber::get_elastance_values(
    std::vector<double>& parameters) {
  double Emax = parameters[global_param_ids[ParamId::EMAX]];
  double Epass = parameters[global_param_ids[ParamId::EPASS]];
  double Vrest = parameters[global_param_ids[ParamId::VREST]];
  double contract_start = parameters[global_param_ids[ParamId::CSTART]];
  double relax_start = parameters[global_param_ids[ParamId::RSTART]];
  double contract_duration = parameters[global_param_ids[ParamId::CDUR]];
  double relax_duration = parameters[global_param_ids[ParamId::RDUR]];

  auto T_HB = model->cardiac_cycle_period;

  double phi = 0;

  auto piecewise_condition = fmod(model->time - contract_start, T_HB);

  if (0 <= piecewise_condition && piecewise_condition < contract_duration) {
    phi = 0.5 * (1 - cos((M_PI * piecewise_condition) / contract_duration));
  } else {
    piecewise_condition = fmod(model->time - relax_start, T_HB);
    if (0 <= piecewise_condition && piecewise_condition < relax_duration) {
      phi = 0.5 * (1 + cos((M_PI * piecewise_condition) / relax_duration));
    }
  }

  Elas = Epass + Emax * phi;
}
