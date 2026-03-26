// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "OpenLoopCoronaryDetailedBC.h"

#include "Model.h"

void OpenLoopCoronaryDetailedBC::setup_dofs(DOFHandler& dofhandler) {
  // 5 equations, 4 internal variables: V_a, V_c, P_a, P_c
  // State variables: P_in, Q_in (from inlet node), V_a, V_c, P_a, P_c
  Block::setup_dofs_(dofhandler, 5, {"V_a", "V_c", "P_a", "P_c"});
}

void OpenLoopCoronaryDetailedBC::update_constant(
    SparseSystem& system, std::vector<double>& parameters) {
  // Parameters
  auto Ra1 = parameters[global_param_ids[0]];
  auto Rv1 = parameters[global_param_ids[1]];
  auto Ca = parameters[global_param_ids[2]];
  auto Cc = parameters[global_param_ids[3]];
  auto Pv = parameters[global_param_ids[5]];
  auto Ra2 = parameters[global_param_ids[6]];

  // Variable indices:
  // 0: P_in, 1: Q_in (from inlet node)
  // 2: V_a, 3: V_c, 4: P_a, 5: P_c (internal variables)

  if (steady) {
    // Steady state equations:
    // 1. V_a = 0
    // 2. V_c = 0
    // 3. P_in - Ra1*Q_in - P_a = 0
    // 4. P_a - (Ra2 + Rv1)*Q_in - P_v = 0 (total resistance from P_a to P_v)
    // 5. P_c - P_v - Rv1*Q_in = 0 (Ohm's law for Rv1)

    system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 1.0;  // V_a = 0
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = 1.0;  // V_c = 0

    // P_in - Ra1*Q_in - P_a = 0
    system.F.coeffRef(global_eqn_ids[2], global_var_ids[0]) = 1.0;
    system.F.coeffRef(global_eqn_ids[2], global_var_ids[1]) = -Ra1;
    system.F.coeffRef(global_eqn_ids[2], global_var_ids[4]) = -1.0;

    // P_a - (Ra2 + Rv1)*Q_in - P_v = 0
    system.F.coeffRef(global_eqn_ids[3], global_var_ids[4]) = 1.0;
    system.F.coeffRef(global_eqn_ids[3], global_var_ids[1]) = -(Ra2 + Rv1);
    system.C(global_eqn_ids[3]) = -Pv;

    // P_c - P_v - Rv1*Q_in = 0
    system.F.coeffRef(global_eqn_ids[4], global_var_ids[5]) = 1.0;
    system.F.coeffRef(global_eqn_ids[4], global_var_ids[1]) = -Rv1;
    system.C(global_eqn_ids[4]) = -Pv;
  } else {
    // Dynamic equations:
    //
    // Eq 0: dV_a/dt = Q_in - (P_a - P_c)/Ra2
    //       dV_a/dt - Q_in + P_a/Ra2 - P_c/Ra2 = 0
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 1.0;     // dV_a/dt
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -1.0;    // -Q_in
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) = 1.0 / Ra2;   // P_a/Ra2
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[5]) = -1.0 / Ra2;  // -P_c/Ra2

    // Eq 1: dV_c/dt = (P_a - P_c)/Ra2 - (P_c - P_v)/Rv1
    //       dV_c/dt - P_a/Ra2 + P_c/Ra2 + P_c/Rv1 - P_v/Rv1 = 0
    system.E.coeffRef(global_eqn_ids[1], global_var_ids[3]) = 1.0;  // dV_c/dt
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[4]) = -1.0 / Ra2;  // -P_a/Ra2
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[5]) =
        1.0 / Ra2 + 1.0 / Rv1;  // P_c*(1/Ra2 + 1/Rv1)

    // Eq 2: P_in - Ra1*Q_in - P_a = 0 (Ohm's law for Ra1)
    system.F.coeffRef(global_eqn_ids[2], global_var_ids[0]) = 1.0;   // P_in
    system.F.coeffRef(global_eqn_ids[2], global_var_ids[1]) = -Ra1;  // -Ra1*Q_in
    system.F.coeffRef(global_eqn_ids[2], global_var_ids[4]) = -1.0;  // -P_a

    // Eq 3: V_a - Ca*P_a = 0 (capacitor relationship for Ca)
    system.F.coeffRef(global_eqn_ids[3], global_var_ids[2]) = 1.0;  // V_a
    system.F.coeffRef(global_eqn_ids[3], global_var_ids[4]) = -Ca;  // -Ca*P_a

    // Eq 4: V_c - Cc*(P_c - Pim) = 0 (capacitor relationship for Cc)
    //       V_c - Cc*P_c + Cc*Pim = 0
    system.F.coeffRef(global_eqn_ids[4], global_var_ids[3]) = 1.0;  // V_c
    system.F.coeffRef(global_eqn_ids[4], global_var_ids[5]) = -Cc;  // -Cc*P_c
    // Note: Cc*Pim is time-dependent and handled in update_time
  }
}

void OpenLoopCoronaryDetailedBC::update_time(SparseSystem& system,
                                              std::vector<double>& parameters) {
  auto Rv1 = parameters[global_param_ids[1]];
  auto Cc = parameters[global_param_ids[3]];
  auto Pim = parameters[global_param_ids[4]];
  auto Pv = parameters[global_param_ids[5]];

  if (!steady) {
    // Update RHS for Eq 1: -P_v/Rv1
    system.C(global_eqn_ids[1]) = -Pv / Rv1;

    // Update RHS for Eq 4: +Cc*Pim
    system.C(global_eqn_ids[4]) = Cc * Pim;
  }
}

void OpenLoopCoronaryDetailedBC::setup_initial_state_dependent_params(
    State initial_state, std::vector<double>& parameters) {
  this->Pim_0 = parameters[global_param_ids[4]];
}
