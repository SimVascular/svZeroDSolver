// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "BloodVesselFC.h"

#include <cmath>

#include "Model.h"

void BloodVesselFC::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {});
}

void BloodVesselFC::update_constant(SparseSystem& system,
                                    std::vector<double>& parameters) {
  // Get parameters
  double inductance = parameters[global_param_ids[ParamId::INDUCTANCE]];
  double resistance = parameters[global_param_ids[ParamId::RESISTANCE]];

  // Get fixed capacitance from model by looking up this block's name
  double capacitance = 0.0;
  std::string block_name = "";
  for (size_t i = 0; i < model->get_num_blocks(); i++) {
    if (model->get_block(i) == this) {
      block_name = model->get_block_name(i);
      break;
    }
  }
  if (!block_name.empty() && model->fixed_capacitance.count(block_name)) {
    capacitance = model->fixed_capacitance.at(block_name);
  }

  // Set element contributions
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[3]) = -inductance;
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[0]) = -capacitance;
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[1]) =
      capacitance * resistance;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -resistance;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
}

void BloodVesselFC::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  // Get fixed capacitance from model
  std::string block_name = "";
  for (size_t i = 0; i < model->get_num_blocks(); i++) {
    if (model->get_block(i) == this) {
      block_name = model->get_block_name(i);
      break;
    }
  }
  double capacitance = 0.0;
  if (!block_name.empty() && model->fixed_capacitance.count(block_name)) {
    capacitance = model->fixed_capacitance.at(block_name);
  }

  double stenosis_coeff = 0.0;
  if (global_param_ids.size() >
      static_cast<size_t>(ParamId::STENOSIS_COEFFICIENT)) {
    stenosis_coeff =
        parameters[global_param_ids[ParamId::STENOSIS_COEFFICIENT]];
  }
  double q_in = y[global_var_ids[1]];
  double dq_in = dy[global_var_ids[1]];
  double stenosis_resistance = stenosis_coeff * fabs(q_in);

  // Set element contributions
  system.C(global_eqn_ids[0]) = stenosis_resistance * -q_in;
  system.C(global_eqn_ids[1]) = stenosis_resistance * 2.0 * capacitance * dq_in;

  double sgn_q_in = (0.0 < q_in) - (q_in < 0.0);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
      stenosis_coeff * sgn_q_in * -2.0 * q_in;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[1]) =
      stenosis_coeff * sgn_q_in * 2.0 * capacitance * dq_in;

  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[1]) =
      stenosis_resistance * 2.0 * capacitance;
}

void BloodVesselFC::update_gradient(
    Eigen::SparseMatrix<double>& jacobian,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& residual,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha, std::vector<double>& y,
    std::vector<double>& dy) {
  // Check if all required observations are available (not NaN)
  auto y0 = y[global_var_ids[0]];
  auto y1 = y[global_var_ids[1]];
  auto y2 = y[global_var_ids[2]];
  auto y3 = y[global_var_ids[3]];

  auto dy0 = dy[global_var_ids[0]];
  auto dy1 = dy[global_var_ids[1]];
  auto dy3 = dy[global_var_ids[3]];

  // Skip residual computation if any required observation is missing (NaN)
  if (std::isnan(y0) || std::isnan(y1) || std::isnan(y2) || std::isnan(y3) ||
      std::isnan(dy0) || std::isnan(dy1) || std::isnan(dy3)) {
    return;
  }

  auto resistance = alpha[global_param_ids[ParamId::RESISTANCE]];
  auto inductance = alpha[global_param_ids[ParamId::INDUCTANCE]];
  double stenosis_coeff = 0.0;

  if (global_param_ids.size() > 2) {
    stenosis_coeff = alpha[global_param_ids[ParamId::STENOSIS_COEFFICIENT]];
  }

  // Get fixed capacitance from model
  std::string block_name = "";
  for (size_t i = 0; i < model->get_num_blocks(); i++) {
    if (model->get_block(i) == this) {
      block_name = model->get_block_name(i);
      break;
    }
  }
  double capacitance = 0.0;
  if (!block_name.empty() && model->fixed_capacitance.count(block_name)) {
    capacitance = model->fixed_capacitance.at(block_name);
  }

  auto stenosis_resistance = stenosis_coeff * fabs(y1);

  // Jacobian entries for R (param 0) - same as BloodVessel
  jacobian.coeffRef(global_eqn_ids[0], global_param_ids[0]) = -y1;
  jacobian.coeffRef(global_eqn_ids[1], global_param_ids[0]) = capacitance * dy1;

  // Jacobian entries for L (param 1) - same as BloodVessel param 2
  jacobian.coeffRef(global_eqn_ids[0], global_param_ids[1]) = -dy3;

  // Note: No Jacobian entries for C since it's fixed!

  if (global_param_ids.size() > 2) {
    // Jacobian entries for stenosis (param 2) - same as BloodVessel param 3
    jacobian.coeffRef(global_eqn_ids[0], global_param_ids[2]) = -fabs(y1) * y1;
    jacobian.coeffRef(global_eqn_ids[1], global_param_ids[2]) =
        2.0 * capacitance * fabs(y1) * dy1;
  }

  // Residual - same as BloodVessel but using fixed capacitance
  residual(global_eqn_ids[0]) =
      y0 - (resistance + stenosis_resistance) * y1 - y2 - inductance * dy3;
  residual(global_eqn_ids[1]) =
      (y1 - y3 - capacitance * dy0 +
       capacitance * (resistance + 2.0 * stenosis_resistance) * dy1);
}
