// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "LinearElastanceChamber.h"

void LinearElastanceChamber::setup_dofs(DOFHandler& dofhandler) {
  // Internal variable is chamber volume
  Block::setup_dofs_(dofhandler, 3, {"Vc"});
}

void LinearElastanceChamber::update_constant(SparseSystem& system,
                                             std::vector<double>& parameters) {
  // Initialize activation function on first call
  if (!activation_func_) {
    initialize_activation_function(parameters);
  }

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

  auto T_cardiac = model->cardiac_cycle_period;
  
  // Compute activation using the activation function
  double phi = activation_func_->compute(model->time, T_cardiac);

  Elas = Epass + Emax * phi;
}

void LinearElastanceChamber::initialize_activation_function(
    std::vector<double>& parameters) {
  // Check if activation_type parameter is provided (optional parameter)
  // Default to PIECEWISE_COSINE (1) for backward compatibility
  int activation_type_int = 1;
  if (global_param_ids.count(ParamId::ACTIVATION_TYPE) > 0) {
    activation_type_int = static_cast<int>(
        parameters[global_param_ids[ParamId::ACTIVATION_TYPE]]);
    activation_params_.type = static_cast<ActivationType>(activation_type_int);
  }

  // Set cardiac period
  activation_params_.cardiac_period = model->cardiac_cycle_period;
  
  // For backward compatibility, if activation parameters weren't set from JSON,
  // try to get them from the parameter vector (for old flat format)
  if (activation_params_.type == ActivationType::PIECEWISE_COSINE) {
    if (global_param_ids.count(ParamId::CSTART) > 0) {
      activation_params_.contract_start = parameters[global_param_ids[ParamId::CSTART]];
    }
    if (global_param_ids.count(ParamId::RSTART) > 0) {
      activation_params_.relax_start = parameters[global_param_ids[ParamId::RSTART]];
    }
    if (global_param_ids.count(ParamId::CDUR) > 0) {
      activation_params_.contract_duration = parameters[global_param_ids[ParamId::CDUR]];
    }
    if (global_param_ids.count(ParamId::RDUR) > 0) {
      activation_params_.relax_duration = parameters[global_param_ids[ParamId::RDUR]];
    }
  }
  
  // Use factory method to create activation function
  activation_func_ = ActivationFunction::create(activation_params_);
}
