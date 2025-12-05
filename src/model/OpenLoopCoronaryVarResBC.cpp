// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "OpenLoopCoronaryVarResBC.h"

#include <cmath>

#include "Model.h"

void OpenLoopCoronaryVarResBC::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {"volume_im"});
}

double OpenLoopCoronaryVarResBC::compute_Ram(double t_cycle, double Ram_min,
                                             double Ram_max, double T_vc,
                                             double T_vr) const {
  // Compute e(t) based on phase in cardiac cycle
  double e_t;
  if (t_cycle <= T_vc) {
    // Contraction phase
    e_t = 0.5 * (1.0 - cos(M_PI * t_cycle / T_vc));
  } else if (t_cycle <= T_vc + T_vr) {
    // Relaxation phase
    e_t = 0.5 * (1.0 + cos(M_PI * (t_cycle - T_vc) / T_vr));
  } else {
    // Rest phase
    e_t = 0.0;
  }

  // Compute time-varying resistance
  double sqrt_Ram_min = sqrt(Ram_min);
  double sqrt_Ram_max = sqrt(Ram_max);
  double Ram_t = sqrt_Ram_min + (sqrt_Ram_max - sqrt_Ram_min) * e_t;
  return Ram_t * Ram_t;
}

void OpenLoopCoronaryVarResBC::update_constant(
    SparseSystem& system, std::vector<double>& parameters) {
  auto Ra = parameters[global_param_ids[0]];
  auto Rv = parameters[global_param_ids[3]];
  auto Ca = parameters[global_param_ids[4]];
  auto Cim = parameters[global_param_ids[5]];

  if (steady) {
    // -P_in + (Ra+Ram+Rv)Q_in + Pv = 0
    // V_im = 0
    // Ram = Ram_min at steady state
    auto Ram_min = parameters[global_param_ids[1]];

    system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = -1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = Ra + Ram_min + Rv;
  } else {
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = Cim * Rv;
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = Cim * Rv;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = -Cim * Rv * Ra;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -Ca * Cim * Rv;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
        Ra * Ca * Cim * Rv;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -Cim * Rv;
  }
}

void OpenLoopCoronaryVarResBC::update_time(SparseSystem& system,
                                           std::vector<double>& parameters) {
  auto Ram_min = parameters[global_param_ids[1]];
  auto Ram_max = parameters[global_param_ids[2]];
  auto Rv = parameters[global_param_ids[3]];
  auto Cim = parameters[global_param_ids[5]];
  auto Pim = parameters[global_param_ids[6]];
  auto Pv = parameters[global_param_ids[7]];
  auto T_vc = parameters[global_param_ids[8]];
  auto T_vr = parameters[global_param_ids[9]];

  // Compute time-varying resistance using cycle time
  double Ram = compute_Ram(model->time, Ram_min, Ram_max, T_vc, T_vr);

  if (steady) {
    system.C(global_eqn_ids[1]) = Pv;
  } else {
    // Update matrix coefficients that depend on Ram
    system.E.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -Cim * Rv * Ram;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -(Rv + Ram);

    // Update C vector
    system.C(global_eqn_ids[0]) =
        Cim * (-Pim + Pv + this->Pim_0 - this->P_Cim_0);
    system.C(global_eqn_ids[1]) =
        (Ram * Cim * Pv) -
        Cim * (Rv + Ram) * (Pim + this->P_Cim_0 - this->Pim_0);
  }
}

void OpenLoopCoronaryVarResBC::setup_initial_state_dependent_params(
    State initial_state, std::vector<double>& parameters) {
  auto P_in = initial_state.y[global_var_ids[0]];
  auto Q_in = initial_state.y[global_var_ids[1]];
  auto P_in_dot = initial_state.ydot[global_var_ids[0]];
  auto Q_in_dot = initial_state.ydot[global_var_ids[1]];
  auto Ra = parameters[global_param_ids[0]];
  auto Ca = parameters[global_param_ids[4]];

  // Ram = Ram_min at initial state
  auto Ram_min = parameters[global_param_ids[1]];

  // Pressure proximal to Ca and distal to Ra
  auto P_Ca = P_in - Ra * Q_in;
  auto P_Ca_dot = P_in_dot - Ra * Q_in_dot;
  // Flow into Ram (inflow minus flow into Ca)
  auto Q_am = Q_in - Ca * P_Ca_dot;
  // Pressure proximal to Cim/Vim and distal to Ram
  this->P_Cim_0 = P_Ca - Ram_min * Q_am;
  // Initial intramyocardial pressure
  this->Pim_0 = parameters[global_param_ids[6]];
}
