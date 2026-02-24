// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberSphere_expmat.h"

#include "Model.h"

void ChamberSphere_expmat::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 7,
                     {"radius", "velo", "stress", "tau", "volume"});
}

void ChamberSphere_expmat::update_constant(SparseSystem& system,
                                    std::vector<double>& parameters) {
  const double rho = parameters[global_param_ids[ParamId::rho]];
  const double thick0 = parameters[global_param_ids[ParamId::thick0]];

  // balance of linear momentum
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[5]) = rho*thick0;

  // spherical stress
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[6]) = -1;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[7]) = 1;

  // volume change
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[8]) = -1;

  // active stress
  system.E.coeffRef(global_eqn_ids[3], global_var_ids[7]) = 1;

  // acceleration
  system.E.coeffRef(global_eqn_ids[4], global_var_ids[4]) = 1;
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[5]) = -1;

  // conservation of mass
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[1]) = 1;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[3]) = -1;
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[8]) = -1;

  // pressure equality
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[0]) = 1;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[2]) = -1;
}

void ChamberSphere_expmat::update_time(SparseSystem& system,
                                std::vector<double>& parameters) {
  // active stress
  get_active_stress_values(parameters);
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[7]) = act;
}

void ChamberSphere_expmat::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  const double sigma_max = parameters[global_param_ids[ParamId::sigma_max]];
  const double C1 = parameters[global_param_ids[ParamId::C1]];
  const double eta = parameters[global_param_ids[ParamId::eta]];
  const double C2 = parameters[global_param_ids[ParamId::C2]];
  const double radius0 = parameters[global_param_ids[ParamId::radius0]];
  const double C0 = parameters[global_param_ids[ParamId::C0]];
  const double C3 = parameters[global_param_ids[ParamId::C3]];
  const double thick0 = parameters[global_param_ids[ParamId::thick0]];

  const double dradius_dt = dy[global_var_ids[4]];
  const double radius = y[global_var_ids[4]];
  const double stress = y[global_var_ids[6]];
  const double Pout = y[global_var_ids[2]];
  const double velo = y[global_var_ids[5]];

  // balance of momentum
  system.C.coeffRef(global_eqn_ids[0]) = (radius + radius0)*(-Pout*(radius + radius0) + stress*thick0)/pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -pow(radius + radius0, 2)/pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) = (-2*Pout*(radius + radius0) + stress*thick0)/pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[6]) = thick0*(radius + radius0)/pow(radius0, 2);

  // spherical stress
  system.C.coeffRef(global_eqn_ids[1]) = 2*(4*C0*C1*(radius + radius0)*(-pow(radius0, 6) + pow(radius + radius0, 6))*(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6))*exp(C1*pow(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6), 2)/(pow(radius0, 4)*pow(radius + radius0, 8))) + 2*C2*C3*pow(radius + radius0, 11)*(-pow(radius0, 2) + pow(radius + radius0, 2))*exp(C3*pow(-pow(radius0, 2) + pow(radius + radius0, 2), 2)/pow(radius0, 4)) + dradius_dt*eta*(2*pow(radius0, 12) + pow(radius + radius0, 12)))/(pow(radius0, 2)*pow(radius + radius0, 11));
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4]) = 2*(32*C0*pow(C1, 2)*(pow(radius0, 6) - pow(radius + radius0, 6))*pow(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6), 2)*(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6) + 3*pow(radius + radius0, 4)*(pow(radius0, 2) - pow(radius + radius0, 2)))*exp(C1*pow(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6), 2)/(pow(radius0, 4)*pow(radius + radius0, 8))) + 8*C2*pow(C3, 2)*pow(radius + radius0, 20)*pow(pow(radius0, 2) - pow(radius + radius0, 2), 2)*exp(C3*pow(pow(radius0, 2) - pow(radius + radius0, 2), 2)/pow(radius0, 4)) + 2*pow(radius0, 4)*pow(radius + radius0, 8)*(12*C0*C1*pow(radius + radius0, 6)*(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6))*exp(C1*pow(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6), 2)/(pow(radius0, 4)*pow(radius + radius0, 8))) + 24*C0*C1*pow(radius + radius0, 4)*(pow(radius0, 2) - pow(radius + radius0, 2))*(pow(radius0, 6) - pow(radius + radius0, 6))*exp(C1*pow(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6), 2)/(pow(radius0, 4)*pow(radius + radius0, 8))) - 2*C0*C1*(pow(radius0, 6) - pow(radius + radius0, 6))*(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6))*exp(C1*pow(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6), 2)/(pow(radius0, 4)*pow(radius + radius0, 8))) + 2*C2*C3*pow(radius + radius0, 12)*exp(C3*pow(pow(radius0, 2) - pow(radius + radius0, 2), 2)/pow(radius0, 4)) - 11*C2*C3*pow(radius + radius0, 10)*(pow(radius0, 2) - pow(radius + radius0, 2))*exp(C3*pow(pow(radius0, 2) - pow(radius + radius0, 2), 2)/pow(radius0, 4)) + 6*dradius_dt*eta*pow(radius + radius0, 11)) + 11*pow(radius0, 4)*pow(radius + radius0, 7)*(4*C0*C1*(radius + radius0)*(pow(radius0, 6) - pow(radius + radius0, 6))*(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6))*exp(C1*pow(pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4) + 2*pow(radius + radius0, 6), 2)/(pow(radius0, 4)*pow(radius + radius0, 8))) + 2*C2*C3*pow(radius + radius0, 11)*(pow(radius0, 2) - pow(radius + radius0, 2))*exp(C3*pow(pow(radius0, 2) - pow(radius + radius0, 2), 2)/pow(radius0, 4)) - dradius_dt*eta*(2*pow(radius0, 12) + pow(radius + radius0, 12))))/(pow(radius0, 6)*pow(radius + radius0, 19));
  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[4]) = 2*eta*(2*pow(radius0, 12) + pow(radius + radius0, 12))/(pow(radius0, 2)*pow(radius + radius0, 11));

  // volume change
  system.C.coeffRef(global_eqn_ids[2]) = 4*M_PI*velo*pow(radius + radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[4]) = 8*M_PI*velo*(radius + radius0);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[5]) = 4*M_PI*pow(radius + radius0, 2);

  // active stress
  system.C.coeffRef(global_eqn_ids[3]) = -act_plus*sigma_max;
}

void ChamberSphere_expmat::get_active_stress_values(std::vector<double>& parameters) {
  const double alpha_max = parameters[global_param_ids[ParamId::alpha_max]];
  const double alpha_min = parameters[global_param_ids[ParamId::alpha_min]];
  const double tsys = parameters[global_param_ids[ParamId::tsys]];
  const double tdias = parameters[global_param_ids[ParamId::tdias]];
  const double steepness = parameters[global_param_ids[ParamId::steepness]];

  const double t = model->time;

  const auto T_cardiac = model->cardiac_cycle_period;
  const auto t_in_cycle = fmod(model->time, T_cardiac);

  const double S_plus = 0.5 * (1.0 + tanh((t_in_cycle - tsys) / steepness));
  const double S_minus = 0.5 * (1.0 - tanh((t_in_cycle - tdias) / steepness));

  // indicator function
  const double f = S_plus * S_minus;

  // activation rates
  const double act_t = alpha_max * f + alpha_min * (1 - f);

  act = std::abs(act_t);
  act_plus = std::max(act_t, 0.0);
}