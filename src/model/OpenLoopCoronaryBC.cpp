// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "OpenLoopCoronaryBC.h"

void OpenLoopCoronaryBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {"volume_im"});
}

void OpenLoopCoronaryBC::update_constant(SparseSystem &system,
                                         std::vector<double> &parameters) {
  auto Ra = parameters[global_param_ids[0]];
  auto Ram = parameters[global_param_ids[1]];
  auto Rv = parameters[global_param_ids[2]];
  auto Ca = parameters[global_param_ids[3]];
  auto Cim = parameters[global_param_ids[4]];

  if (steady) {
    // Different equations for steady initial condition
    // Equations are:
    // -P_in + (Ra+Ram+Rv)Q_in + Pv = 0
    // V_im = 0
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = -1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = Ra + Ram + Rv;
  } else {
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = Cim * Rv;
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = Cim * Rv;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = -Cim * Rv * Ra;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -(Rv + Ram);

    system.E.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -Ca * Cim * Rv;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
        Ra * Ca * Cim * Rv;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -Cim * Rv;
    system.E.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -Cim * Rv * Ram;
  }
}

void OpenLoopCoronaryBC::update_time(SparseSystem &system,
                                     std::vector<double> &parameters) {
  auto Ram = parameters[global_param_ids[1]];
  auto Rv = parameters[global_param_ids[2]];
  auto Cim = parameters[global_param_ids[4]];
  auto Pim = parameters[global_param_ids[5]];
  auto Pv = parameters[global_param_ids[6]];

  if (steady) {
    system.C(global_eqn_ids[1]) = Pv;
  } else {
    system.C(global_eqn_ids[0]) =
        Cim * (-Pim + Pv + this->Pim_0 - this->P_Cim_0);
    system.C(global_eqn_ids[1]) =
        (Ram * Cim * Pv) -
        Cim * (Rv + Ram) * (Pim + this->P_Cim_0 - this->Pim_0);
  }
}

void OpenLoopCoronaryBC::setup_initial_state_dependent_params(
    State initial_state, std::vector<double> &parameters) {
  auto P_in = initial_state.y[global_var_ids[0]];
  auto Q_in = initial_state.y[global_var_ids[1]];
  auto P_in_dot = initial_state.ydot[global_var_ids[0]];
  auto Q_in_dot = initial_state.ydot[global_var_ids[1]];
  auto Ra = parameters[global_param_ids[0]];
  auto Ram = parameters[global_param_ids[1]];
  auto Ca = parameters[global_param_ids[3]];
  // Pressure proximal to Ca and distal to Ra
  auto P_Ca = P_in - Ra * Q_in;
  auto P_Ca_dot = P_in_dot - Ra * Q_in_dot;
  // Flow into Ram (inflow minus flow into Ca)
  auto Q_am = Q_in - Ca * P_Ca_dot;
  // Pressure proximal to Cim/Vim and distal to Ram
  this->P_Cim_0 = P_Ca - Ram * Q_am;
  // Initial intramyocardial pressure
  this->Pim_0 = parameters[global_param_ids[5]];
}
