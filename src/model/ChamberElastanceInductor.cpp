// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ChamberElastanceInductor.h"

void ChamberElastanceInductor::setup_dofs(DOFHandler& dofhandler) {
  // Internal variable is chamber volume
  Block::setup_dofs_(dofhandler, 3, {"Vc"});
}

void ChamberElastanceInductor::update_constant(
    SparseSystem& system, std::vector<double>& parameters) {
  // Initialize activation function on first call
  if (!activation_func_) {
    initialize_activation_function(parameters);
  }

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

  // Eq 0: P_in - E(t)(Vc - Vrest) = P_in - E(t)*Vc + E(t)*Vrest = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -1 * Elas;
  system.C.coeffRef(global_eqn_ids[0]) = Elas * Vrest;
}

void ChamberElastanceInductor::get_elastance_values(
    std::vector<double>& parameters) {
  double Emax = parameters[global_param_ids[ParamId::EMAX]];
  double Emin = parameters[global_param_ids[ParamId::EMIN]];
  double Vrd = parameters[global_param_ids[ParamId::VRD]];
  double Vrs = parameters[global_param_ids[ParamId::VRS]];

  auto T_cardiac = model->cardiac_cycle_period;
  
  // Compute activation using the activation function
  double act = activation_func_->compute(model->time, T_cardiac);

  Vrest = (1.0 - act) * (Vrd - Vrs) + Vrs;
  Elas = (Emax - Emin) * act + Emin;
}

void ChamberElastanceInductor::initialize_activation_function(
    std::vector<double>& parameters) {
  // Check if activation_type parameter is provided (optional parameter)
  int activation_type_int = 0;  // Default to HALF_COSINE
  if (global_param_ids.count(ParamId::ACTIVATION_TYPE) > 0) {
    activation_type_int = static_cast<int>(
        parameters[global_param_ids[ParamId::ACTIVATION_TYPE]]);
  }

  auto T_cardiac = model->cardiac_cycle_period;

  switch (activation_type_int) {
    case 0:  // HALF_COSINE (default)
    {
      double t_active = parameters[global_param_ids[ParamId::TACTIVE]];
      double t_twitch = parameters[global_param_ids[ParamId::TTWITCH]];
      activation_func_ = std::make_unique<HalfCosineActivation>(t_active, t_twitch);
      break;
    }
    case 1:  // PIECEWISE_COSINE
    {
      double contract_start = parameters[global_param_ids[ParamId::CSTART]];
      double relax_start = parameters[global_param_ids[ParamId::RSTART]];
      double contract_duration = parameters[global_param_ids[ParamId::CDUR]];
      double relax_duration = parameters[global_param_ids[ParamId::RDUR]];
      activation_func_ = std::make_unique<PiecewiseCosineActivation>(
          contract_start, relax_start, contract_duration, relax_duration);
      break;
    }
    case 2:  // TWO_HILL
    {
      double t_shift = parameters[global_param_ids[ParamId::TSHIFT]];
      double tau_1 = parameters[global_param_ids[ParamId::TAU_1]];
      double tau_2 = parameters[global_param_ids[ParamId::TAU_2]];
      double m1 = parameters[global_param_ids[ParamId::M1]];
      double m2 = parameters[global_param_ids[ParamId::M2]];
      activation_func_ = std::make_unique<TwoHillActivation>(
          t_shift, tau_1, tau_2, m1, m2, T_cardiac);
      break;
    }
    default:
      throw std::runtime_error(
          "ChamberElastanceInductor: Invalid activation_type " +
          std::to_string(activation_type_int));
  }
}
