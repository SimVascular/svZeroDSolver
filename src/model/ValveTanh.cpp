// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ValveTanh.h"

#include "Model.h"

void ValveTanh::setup_dofs(DOFHandler& dofhandler) {
  // 2 eqns, one for pressure and one for flow. The valve status feeds back
  // into nothing, so it is a post-processing variable rather than an unknown
  Block::setup_dofs_(dofhandler, 2, {}, {"valve_status"});
}

// update_constant updates matrices E and F from E(y,t)*y_dot + F(y,t)*y +
// c(y,t) = 0 with terms that DO NOT DEPEND ON THE SOLUTION
void ValveTanh::update_constant(SparseSystem& system,
                                std::vector<double>& parameters) {
  // Set element contributions
  // coeffRef args are the indices (i,j) of the matrix
  // global_eqn_ids: number of rows in the matrix, set in setup_dofs
  // global_var_ids: number of columns, organized as pressure and flow of all
  // inlets and then all outlets of the block
  double Rmin = parameters[global_param_ids[ParamId::RMIN]];
  double Rmax = parameters[global_param_ids[ParamId::RMAX]];
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
      -0.5 * (Rmax + Rmin);
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
}

void ValveTanh::update_outputs(
    Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& ydot) {
  double p_in = y[global_var_ids[0]];
  double p_out = y[global_var_ids[2]];
  double steep =
      model->get_parameter_value(global_param_ids[ParamId::STEEPNESS]);
  y[global_var_ids[4]] = 0.5 * (1.0 + tanh(steep * (p_out - p_in)));
}

// update_solution updates matrices E and F from E(y,t)*y_dot + F(y,t)*y +
// c(y,t) = 0 with terms that DO DEPEND ON THE SOLUTION (will change with each
// time step)
void ValveTanh::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  // Get states
  double p_in = y[global_var_ids[0]];
  double p_out = y[global_var_ids[2]];
  double q_in = y[global_var_ids[1]];
  // Get parameters
  double Rmin = parameters[global_param_ids[ParamId::RMIN]];
  double Rmax = parameters[global_param_ids[ParamId::RMAX]];
  double steep = parameters[global_param_ids[ParamId::STEEPNESS]];

  // Helper function
  double fun_tanh = tanh(steep * (p_out - p_in));

  // Nonlinear term
  system.C(global_eqn_ids[0]) = -0.5 * q_in * (Rmax - Rmin) * fun_tanh;

  // Derivatives of non-linear terms
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[0]) =
      0.5 * q_in * (Rmax - Rmin) * steep * (1.0 - pow(fun_tanh, 2));
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
      -0.5 * (Rmax - Rmin) * fun_tanh;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) =
      -0.5 * q_in * (Rmax - Rmin) * steep * (1.0 - pow(fun_tanh, 2));
}
