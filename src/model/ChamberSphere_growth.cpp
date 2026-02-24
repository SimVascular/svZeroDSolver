// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberSphere_growth.h"

#include "Model.h"

void ChamberSphere_growth::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 7,
                     {"radius", "velo", "stress", "tau", "volume"});
}

void ChamberSphere_growth::update_constant(SparseSystem& system,
                                    std::vector<double>& parameters) {
  const double Theta_r = parameters[global_param_ids[ParamId::Theta_r]];
  const double rho = parameters[global_param_ids[ParamId::rho]];
  const double Theta_c = parameters[global_param_ids[ParamId::Theta_c]];
  const double thick0 = parameters[global_param_ids[ParamId::thick0]];
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[5]) = pow(Theta_c, 2)*Theta_r*rho*thick0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[6]) = -1;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[7]) = 1;
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[8]) = -1;
  system.E.coeffRef(global_eqn_ids[3], global_var_ids[7]) = 1;
  system.E.coeffRef(global_eqn_ids[4], global_var_ids[4]) = 1;
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[5]) = -1;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[1]) = 1;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[3]) = -1;
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[8]) = -1;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[0]) = 1;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[2]) = -1;

}

void ChamberSphere_growth::update_time(SparseSystem& system,
                                std::vector<double>& parameters) {
  // active stress
  compute_active_stress_values(parameters);
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[7]) = act;
}

void ChamberSphere_growth::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
        const double C1 = parameters[global_param_ids[ParamId::C1]];
        const double Theta_r = parameters[global_param_ids[ParamId::Theta_r]];
        const double C3 = parameters[global_param_ids[ParamId::C3]];
        const double sigma_max = parameters[global_param_ids[ParamId::sigma_max]];
        const double C2 = parameters[global_param_ids[ParamId::C2]];
        const double C0 = parameters[global_param_ids[ParamId::C0]];
        const double thick0 = parameters[global_param_ids[ParamId::thick0]];
        const double eta = parameters[global_param_ids[ParamId::eta]];
        const double radius0 = parameters[global_param_ids[ParamId::radius0]];
        const double Theta_c = parameters[global_param_ids[ParamId::Theta_c]];
        // compute time dependent constant act_plus
        const double velo = y[global_var_ids[5]];
        const double dradius_dt = dy[global_var_ids[4]];
        const double stress = y[global_var_ids[6]];
        const double Pout = y[global_var_ids[2]];
        const double radius = y[global_var_ids[4]];
        system.C.coeffRef(global_eqn_ids[0]) = Theta_c*(radius + radius0)*(-Pout*Theta_c*(radius + radius0) + stress*thick0)/pow(radius0, 2);
        system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -pow(Theta_c, 2)*pow(radius + radius0, 2)/pow(radius0, 2);
        system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) = Theta_c*(-2*Pout*Theta_c*(radius + radius0) + stress*thick0)/pow(radius0, 2);
        system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[6]) = Theta_c*thick0*(radius + radius0)/pow(radius0, 2);
        system.C.coeffRef(global_eqn_ids[1]) = 2*(-4*C0*C1*(radius + radius0)*(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)))*(pow(Theta_r, 2)*pow(radius0, 6) - pow(Theta_r, 2)*pow(radius + radius0, 6) + (pow(Theta_c, 2)*pow(Theta_r, 2) - 1)*(0.66666666666666663*pow(Theta_c, 2)*pow(radius + radius0, 6) + 0.33333333333333331*pow(Theta_r, 2)*pow(radius0, 6)))*exp(C1*pow(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(radius + radius0, 8)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) + C2*C3*(1.3333333333333333*pow(Theta_c, 2) + 0.66666666666666663*pow(Theta_r, 2))*pow(radius + radius0, 11)*(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)))*exp(C3*pow(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) + pow(Theta_r, 2)*dradius_dt*eta*(pow(Theta_c, 2)*pow(radius + radius0, 12) + 2*pow(Theta_r, 2)*pow(radius0, 12))*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))/(pow(Theta_r, 2)*pow(radius0, 2)*pow(radius + radius0, 11)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2));
        system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4]) = 2*(-32*C0*pow(C1, 2)*pow(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)*(pow(Theta_r, 2)*pow(radius0, 6) - pow(Theta_r, 2)*pow(radius + radius0, 6) + (pow(Theta_c, 2)*pow(Theta_r, 2) - 1)*(0.66666666666666663*pow(Theta_c, 2)*pow(radius + radius0, 6) + 0.33333333333333331*pow(Theta_r, 2)*pow(radius0, 6)))*(-2*pow(Theta_c, 2)*pow(radius + radius0, 6) - pow(Theta_r, 2)*pow(radius0, 6) + 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)) + 3*pow(radius + radius0, 4)*(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663))))*exp(C1*pow(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(radius + radius0, 8)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) + 4*C2*pow(C3, 2)*pow(Theta_c, 2)*(1.3333333333333333*pow(Theta_c, 2) + 0.66666666666666663*pow(Theta_r, 2))*pow(radius + radius0, 20)*pow(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)*exp(C3*pow(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) + pow(radius0, 4)*pow(radius + radius0, 8)*(-4*C0*C1*pow(radius + radius0, 6)*(4.0*pow(Theta_c, 2)*(pow(Theta_c, 2)*pow(Theta_r, 2) - 1) - 6*pow(Theta_r, 2))*(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)))*exp(C1*pow(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(radius + radius0, 8)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) - 48*C0*C1*pow(radius + radius0, 4)*(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)))*(pow(Theta_r, 2)*pow(radius0, 6) - pow(Theta_r, 2)*pow(radius + radius0, 6) + (pow(Theta_c, 2)*pow(Theta_r, 2) - 1)*(0.66666666666666663*pow(Theta_c, 2)*pow(radius + radius0, 6) + 0.33333333333333331*pow(Theta_r, 2)*pow(radius0, 6)))*exp(C1*pow(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(radius + radius0, 8)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) - 4*C0*C1*(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)))*(pow(Theta_r, 2)*pow(radius0, 6) - pow(Theta_r, 2)*pow(radius + radius0, 6) + (pow(Theta_c, 2)*pow(Theta_r, 2) - 1)*(0.66666666666666663*pow(Theta_c, 2)*pow(radius + radius0, 6) + 0.33333333333333331*pow(Theta_r, 2)*pow(radius0, 6)))*exp(C1*pow(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(radius + radius0, 8)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) + 2*C2*C3*pow(Theta_c, 2)*(1.3333333333333333*pow(Theta_c, 2) + 0.66666666666666663*pow(Theta_r, 2))*pow(radius + radius0, 12)*exp(C3*pow(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) + 11*C2*C3*(1.3333333333333333*pow(Theta_c, 2) + 0.66666666666666663*pow(Theta_r, 2))*pow(radius + radius0, 10)*(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)))*exp(C3*pow(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) + 12*pow(Theta_c, 2)*pow(Theta_r, 2)*dradius_dt*eta*pow(radius + radius0, 11)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2) + 11*pow(radius0, 4)*pow(radius + radius0, 7)*(4*C0*C1*(radius + radius0)*(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)))*(pow(Theta_r, 2)*pow(radius0, 6) - pow(Theta_r, 2)*pow(radius + radius0, 6) + (pow(Theta_c, 2)*pow(Theta_r, 2) - 1)*(0.66666666666666663*pow(Theta_c, 2)*pow(radius + radius0, 6) + 0.33333333333333331*pow(Theta_r, 2)*pow(radius0, 6)))*exp(C1*pow(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(radius + radius0, 8)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) - C2*C3*(1.3333333333333333*pow(Theta_c, 2) + 0.66666666666666663*pow(Theta_r, 2))*pow(radius + radius0, 11)*(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)))*exp(C3*pow(pow(Theta_c, 2)*pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/(pow(radius0, 4)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) - pow(Theta_r, 2)*dradius_dt*eta*(pow(Theta_c, 2)*pow(radius + radius0, 12) + 2*pow(Theta_r, 2)*pow(radius0, 12))*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))/(pow(Theta_r, 2)*pow(radius0, 6)*pow(radius + radius0, 19)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 4));
        system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[4]) = 2*eta*(pow(Theta_c, 2)*pow(radius + radius0, 12) + 2*pow(Theta_r, 2)*pow(radius0, 12))/(pow(radius0, 2)*pow(radius + radius0, 11));
        system.C.coeffRef(global_eqn_ids[2]) = 4*M_PI*velo*pow(radius + radius0, 2);
        system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[4]) = 8*M_PI*velo*(radius + radius0);
        system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[5]) = 4*M_PI*pow(radius + radius0, 2);
        system.C.coeffRef(global_eqn_ids[3]) = -act_plus*sigma_max;
}

void ChamberSphere_growth::compute_active_stress_values(std::vector<double>& parameters) {
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