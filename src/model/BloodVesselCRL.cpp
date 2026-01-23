// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "BloodVesselCRL.h"

void BloodVesselCRL::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {});
}

void BloodVesselCRL::update_constant(SparseSystem& system,
                                     std::vector<double>& parameters) {
  // Get parameters
  double capacitance = parameters[global_param_ids[ParamId::CAPACITANCE]];
  double inductance = parameters[global_param_ids[ParamId::INDUCTANCE]];
  double resistance = parameters[global_param_ids[ParamId::RESISTANCE]];

  // Set element contributions
  // coeffRef args are the indices (i,j) of the matrix
  // global_eqn_ids: number of rows in the matrix, set in setup_dofs
  // global_var_ids: number of columns, organized as pressure and flow of all
  // inlets and then all outlets of the block
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[3]) = -inductance;
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[0]) = -capacitance;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[3]) = -resistance;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
}

void BloodVesselCRL::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  // Get parameters
  double capacitance = parameters[global_param_ids[ParamId::CAPACITANCE]];
  double stenosis_coeff =
      parameters[global_param_ids[ParamId::STENOSIS_COEFFICIENT]];
  double q_out = y[global_var_ids[3]];
  double dq_out = dy[global_var_ids[3]];
  double stenosis_resistance = stenosis_coeff * fabs(q_out);

  // Set element contributions
  system.C(global_eqn_ids[0]) = stenosis_resistance * -q_out;

  double sgn_q_out = (0.0 < q_out) - (q_out < 0.0);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
      stenosis_coeff * sgn_q_out * -2.0 * q_out;
}

void BloodVesselCRL::update_gradient(
    Eigen::SparseMatrix<double>& jacobian,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& residual,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha, std::vector<double>& y,
    std::vector<double>& dy) {
  auto y0 = y[global_var_ids[0]];
  auto y1 = y[global_var_ids[1]];
  auto y2 = y[global_var_ids[2]];
  auto y3 = y[global_var_ids[3]];

  auto dy0 = dy[global_var_ids[0]];
  auto dy1 = dy[global_var_ids[1]];
  auto dy3 = dy[global_var_ids[3]];

  auto resistance = alpha[global_param_ids[ParamId::RESISTANCE]];
  auto capacitance = alpha[global_param_ids[ParamId::CAPACITANCE]];
  auto inductance = alpha[global_param_ids[ParamId::INDUCTANCE]];
  double stenosis_coeff = 0.0;

  if (global_param_ids.size() > 3) {
    stenosis_coeff = alpha[global_param_ids[ParamId::STENOSIS_COEFFICIENT]];
  }
  auto stenosis_resistance = stenosis_coeff * fabs(y3);

  jacobian.coeffRef(global_eqn_ids[0], global_param_ids[0]) = -y3;
  jacobian.coeffRef(global_eqn_ids[0], global_param_ids[2]) = -dy3;

  if (global_param_ids.size() > 3) {
    jacobian.coeffRef(global_eqn_ids[0], global_param_ids[3]) = -fabs(y3) * y3;
  }

  jacobian.coeffRef(global_eqn_ids[1], global_param_ids[1]) = -dy0;

  residual(global_eqn_ids[0]) =
      y0 - (resistance + stenosis_resistance) * y3 - y2 - inductance * dy3;
  residual(global_eqn_ids[1]) = y1 - y3 - capacitance * dy0;
}
