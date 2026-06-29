// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberSphere_StrainDepActStress.h"

#include "Model.h"

void ChamberSphere_StrainDepActStress::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 11,
                     {"radius", "velo", "stress", "tau", "volume", "e_c",
                      "tau_c", "k_c", "omega"});
}

void ChamberSphere_StrainDepActStress::update_constant(
    SparseSystem& system, std::vector<double>& parameters) {
  const double alpha_r = parameters[global_param_ids[ParamId::alpha_r]];
  const double thick0 = parameters[global_param_ids[ParamId::thick0]];
  const double mu = parameters[global_param_ids[ParamId::mu]];
  const double rho = parameters[global_param_ids[ParamId::rho]];

  // balance of linear momentum
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[5]) = rho * thick0;

  // spherical stress
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[6]) = -1;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[7]) = 1;

  // volume change
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[8]) = -1;

  // acceleration
  system.E.coeffRef(global_eqn_ids[3], global_var_ids[4]) = 1;
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[5]) = -1;

  // conservation of mass
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[1]) = 1;
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[3]) = -1;
  system.E.coeffRef(global_eqn_ids[4], global_var_ids[8]) = -1;

  // pressure equality
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[0]) = 1;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[2]) = -1;

  // active stress
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[7]) = -1;

  // relaxation dynamics (d_omega_dt)
  system.E.coeffRef(global_eqn_ids[7], global_var_ids[12]) = -alpha_r;
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[12]) = -1;

  // sarcomere strain (d_e_c_dt)
  system.E.coeffRef(global_eqn_ids[8], global_var_ids[9]) = -mu;
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[10]) = -1;

  // sarcomere stiffness (d_k_c_dt)
  system.E.coeffRef(global_eqn_ids[9], global_var_ids[11]) = -1;

  // sarcomere active stress (d_tau_c_dt)
  system.E.coeffRef(global_eqn_ids[10], global_var_ids[10]) = -1;
}

void ChamberSphere_StrainDepActStress::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  get_active_stress_values(parameters, y[global_var_ids[9]]);

  const double eta = parameters[global_param_ids[ParamId::eta]];
  const double C0 = parameters[global_param_ids[ParamId::C0]];
  const double C1 = parameters[global_param_ids[ParamId::C1]];
  const double C2 = parameters[global_param_ids[ParamId::C2]];
  const double C3 = parameters[global_param_ids[ParamId::C3]];
  const double alpha = parameters[global_param_ids[ParamId::alpha]];
  const double radius0 = parameters[global_param_ids[ParamId::radius0]];
  const double k_0 = parameters[global_param_ids[ParamId::k_0]];
  const double sigma_0 = parameters[global_param_ids[ParamId::sigma_0]];
  const double E_s = parameters[global_param_ids[ParamId::E_s]];
  const double thick0 = parameters[global_param_ids[ParamId::thick0]];

  const double k_c = y[global_var_ids[11]];
  const double velo = y[global_var_ids[5]];
  const double e_c = y[global_var_ids[9]];
  const double Pout = y[global_var_ids[2]];
  const double dradius_dt = dy[global_var_ids[4]];
  const double radius = y[global_var_ids[4]];
  const double stress = y[global_var_ids[6]];
  const double omega = y[global_var_ids[12]];
  const double tau_c = y[global_var_ids[10]];
  const double de_c_dt = dy[global_var_ids[9]];

  // balance of linear momentum
  system.C.coeffRef(global_eqn_ids[0]) =
      (radius + radius0) * (-Pout * (radius + radius0) + stress * thick0) /
      pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) =
      -pow(radius + radius0, 2) / pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      (-2 * Pout * (radius + radius0) + stress * thick0) / pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[6]) =
      thick0 * (radius + radius0) / pow(radius0, 2);

  // spherical stress
  system.C.coeffRef(global_eqn_ids[1]) =
      2 *
      (4 * C0 * C1 * (radius + radius0) *
           (-pow(radius0, 6) + pow(radius + radius0, 6)) *
           (pow(radius0, 6) - 3 * pow(radius0, 2) * pow(radius + radius0, 4) +
            2 * pow(radius + radius0, 6)) *
           exp(C1 *
               pow(pow(radius0, 6) -
                       3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                       2 * pow(radius + radius0, 6),
                   2) /
               (pow(radius0, 4) * pow(radius + radius0, 8))) +
       2 * C2 * C3 * pow(radius + radius0, 11) *
           (-pow(radius0, 2) + pow(radius + radius0, 2)) *
           exp(C3 * pow(-pow(radius0, 2) + pow(radius + radius0, 2), 2) /
               pow(radius0, 4)) +
       dradius_dt * eta * (2 * pow(radius0, 12) + pow(radius + radius0, 12))) /
      (pow(radius0, 2) * pow(radius + radius0, 11));
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4]) =
      2 *
      (32 * C0 * pow(C1, 2) * (pow(radius0, 6) - pow(radius + radius0, 6)) *
           pow(pow(radius0, 6) -
                   3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                   2 * pow(radius + radius0, 6),
               2) *
           (pow(radius0, 6) - 3 * pow(radius0, 2) * pow(radius + radius0, 4) +
            2 * pow(radius + radius0, 6) +
            3 * pow(radius + radius0, 4) *
                (pow(radius0, 2) - pow(radius + radius0, 2))) *
           exp(C1 *
               pow(pow(radius0, 6) -
                       3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                       2 * pow(radius + radius0, 6),
                   2) /
               (pow(radius0, 4) * pow(radius + radius0, 8))) +
       8 * C2 * pow(C3, 2) * pow(radius + radius0, 20) *
           pow(pow(radius0, 2) - pow(radius + radius0, 2), 2) *
           exp(C3 * pow(pow(radius0, 2) - pow(radius + radius0, 2), 2) /
               pow(radius0, 4)) +
       2 * pow(radius0, 4) * pow(radius + radius0, 8) *
           (12 * C0 * C1 * pow(radius + radius0, 6) *
                (pow(radius0, 6) -
                 3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                 2 * pow(radius + radius0, 6)) *
                exp(C1 *
                    pow(pow(radius0, 6) -
                            3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                            2 * pow(radius + radius0, 6),
                        2) /
                    (pow(radius0, 4) * pow(radius + radius0, 8))) +
            24 * C0 * C1 * pow(radius + radius0, 4) *
                (pow(radius0, 2) - pow(radius + radius0, 2)) *
                (pow(radius0, 6) - pow(radius + radius0, 6)) *
                exp(C1 *
                    pow(pow(radius0, 6) -
                            3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                            2 * pow(radius + radius0, 6),
                        2) /
                    (pow(radius0, 4) * pow(radius + radius0, 8))) -
            2 * C0 * C1 * (pow(radius0, 6) - pow(radius + radius0, 6)) *
                (pow(radius0, 6) -
                 3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                 2 * pow(radius + radius0, 6)) *
                exp(C1 *
                    pow(pow(radius0, 6) -
                            3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                            2 * pow(radius + radius0, 6),
                        2) /
                    (pow(radius0, 4) * pow(radius + radius0, 8))) +
            2 * C2 * C3 * pow(radius + radius0, 12) *
                exp(C3 * pow(pow(radius0, 2) - pow(radius + radius0, 2), 2) /
                    pow(radius0, 4)) -
            11 * C2 * C3 * pow(radius + radius0, 10) *
                (pow(radius0, 2) - pow(radius + radius0, 2)) *
                exp(C3 * pow(pow(radius0, 2) - pow(radius + radius0, 2), 2) /
                    pow(radius0, 4)) +
            6 * dradius_dt * eta * pow(radius + radius0, 11)) +
       11 * pow(radius0, 4) * pow(radius + radius0, 7) *
           (4 * C0 * C1 * (radius + radius0) *
                (pow(radius0, 6) - pow(radius + radius0, 6)) *
                (pow(radius0, 6) -
                 3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                 2 * pow(radius + radius0, 6)) *
                exp(C1 *
                    pow(pow(radius0, 6) -
                            3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                            2 * pow(radius + radius0, 6),
                        2) /
                    (pow(radius0, 4) * pow(radius + radius0, 8))) +
            2 * C2 * C3 * pow(radius + radius0, 11) *
                (pow(radius0, 2) - pow(radius + radius0, 2)) *
                exp(C3 * pow(pow(radius0, 2) - pow(radius + radius0, 2), 2) /
                    pow(radius0, 4)) -
            dradius_dt * eta *
                (2 * pow(radius0, 12) + pow(radius + radius0, 12)))) /
      (pow(radius0, 6) * pow(radius + radius0, 19));
  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[4]) =
      2 * eta * (2 * pow(radius0, 12) + pow(radius + radius0, 12)) /
      (pow(radius0, 2) * pow(radius + radius0, 11));

  // volume change
  system.C.coeffRef(global_eqn_ids[2]) =
      4 * M_PI * velo * pow(radius + radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[4]) =
      8 * M_PI * velo * (radius + radius0);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[5]) =
      4 * M_PI * pow(radius + radius0, 2);

  // active stress
  system.C.coeffRef(global_eqn_ids[6]) =
      E_s * (-e_c * pow(radius0, 2) + 0.5 * pow(radius, 2) + radius * radius0) /
      (pow(radius0, 2) * pow(2 * e_c + 1, 2));
  system.dC_dy.coeffRef(global_eqn_ids[6], global_var_ids[4]) =
      E_s * (1.0 * radius + radius0) / (pow(radius0, 2) * pow(2 * e_c + 1, 2));
  system.dC_dy.coeffRef(global_eqn_ids[6], global_var_ids[9]) =
      E_s *
      (4 * e_c * pow(radius0, 2) - 2.0 * pow(radius, 2) - 4 * radius * radius0 -
       pow(radius0, 2) * (2 * e_c + 1)) /
      (pow(radius0, 2) * pow(2 * e_c + 1, 3));

  // relaxation dynamics (d_omega_dt)
  system.C.coeffRef(global_eqn_ids[7]) = m_0;

  // sarcomere element strain (d_e_c_dt)
  system.C.coeffRef(global_eqn_ids[8]) =
      E_s * (pow(radius, 2) + 2 * radius * radius0 + pow(radius0, 2)) *
      (-e_c * pow(radius0, 2) + 0.5 * pow(radius, 2) + radius * radius0) /
      (pow(radius0, 4) * pow(2 * e_c + 1, 3));
  system.dC_dy.coeffRef(global_eqn_ids[8], global_var_ids[4]) =
      E_s *
      (2 * (radius + radius0) *
           (-e_c * pow(radius0, 2) + 0.5 * pow(radius, 2) + radius * radius0) +
       (1.0 * radius + radius0) *
           (pow(radius, 2) + 2 * radius * radius0 + pow(radius0, 2))) /
      (pow(radius0, 4) * pow(2 * e_c + 1, 3));
  system.dC_dy.coeffRef(global_eqn_ids[8], global_var_ids[9]) =
      E_s * (pow(radius, 2) + 2 * radius * radius0 + pow(radius0, 2)) *
      (6 * e_c * pow(radius0, 2) - 3.0 * pow(radius, 2) - 6 * radius * radius0 -
       pow(radius0, 2) * (2 * e_c + 1)) /
      (pow(radius0, 4) * pow(2 * e_c + 1, 4));

  // sarcomere stiffness (d_k_c_dt)
  system.C.coeffRef(global_eqn_ids[9]) =
      k_0 * n_0 * u_plus -
      k_c * (alpha * fabs(de_c_dt) + omega * fabs(u_minus) + u_plus);
  system.dC_dydot.coeffRef(global_eqn_ids[9], global_var_ids[9]) =
      ((de_c_dt == 0) ? (0) : (-alpha * de_c_dt * k_c / fabs(de_c_dt)));
  system.dC_dy.coeffRef(global_eqn_ids[9], global_var_ids[11]) =
      -alpha * fabs(de_c_dt) - omega * fabs(u_minus) - u_plus;
  system.dC_dy.coeffRef(global_eqn_ids[9], global_var_ids[12]) =
      -k_c * fabs(u_minus);

  // sarcomere active stress (d_tau_c_dt)
  system.C.coeffRef(global_eqn_ids[10]) =
      de_c_dt * k_c + n_0 * sigma_0 * u_plus -
      tau_c * (alpha * fabs(de_c_dt) + omega * fabs(u_minus) + u_plus);
  system.dC_dydot.coeffRef(global_eqn_ids[10], global_var_ids[9]) =
      ((de_c_dt == 0) ? (k_c)
                      : (-alpha * de_c_dt * tau_c / fabs(de_c_dt) + k_c));
  system.dC_dy.coeffRef(global_eqn_ids[10], global_var_ids[10]) =
      -alpha * fabs(de_c_dt) - omega * fabs(u_minus) - u_plus;
  system.dC_dy.coeffRef(global_eqn_ids[10], global_var_ids[11]) = de_c_dt;
  system.dC_dy.coeffRef(global_eqn_ids[10], global_var_ids[12]) =
      -tau_c * fabs(u_minus);
}

void ChamberSphere_StrainDepActStress::get_active_stress_values(
    std::vector<double>& parameters, const double e_c) {
  //   const double E_s = parameters[global_param_ids[ParamId::E_s]];
  //   const double mu = parameters[global_param_ids[ParamId::mu]];
  //   const double alpha_r = parameters[global_param_ids[ParamId::alpha_r]];
  //   const double alpha = parameters[global_param_ids[ParamId::alpha]];
  //   const double k_0 = parameters[global_param_ids[ParamId::k_0]];
  //   const double sigma_0 = parameters[global_param_ids[ParamId::sigma_0]];

  const double t = model->time;

  const auto T_cardiac = model->cardiac_cycle_period;
  const auto t_in_cycle = fmod(model->time, T_cardiac);

  n_0 = 0.0;
  m_0 = 0.0;
  double u = 0.0;
  u_plus = 0.0;
  u_minus = 0.0;

  // activation strain-dependence
  if (e_c > -0.4 && e_c <= 0.3) {
    n_0 = 0.22 + 0.53 * e_c;
  } else if (e_c > 0.3 && e_c <= 1.0) {
    n_0 = 0.112857125 + 0.8871428571 * e_c;
  } else if (e_c > 1.0 && e_c <= 1.3) {
    n_0 = 1;
  } else if (e_c > 1.3 && e_c <= 2.4) {
    n_0 = 2.182 - 0.9091 * e_c;
  } else {
    n_0 = 0.0;
  }

  // relaxation strain-dependence
  if (e_c <= 0.95) {
    m_0 = 1.87;
  } else if (e_c > 0.95 && e_c <= 1.0) {
    m_0 = 18.4 - 17.4 * e_c;
  } else {
    m_0 = 1.0;
  }

  // activation reaction rate function
  if (t_in_cycle <= 0.1) {
    u = -3.2;
  } else if (t_in_cycle <= 0.27 && t_in_cycle > 0.1) {
    u = -3.2 + 38.2 / 0.17 * (t_in_cycle-0.1);
  } else if (t_in_cycle > 0.27 && t_in_cycle <= 0.32) {
    u = 35;
  } else if (t_in_cycle > 0.32 && t_in_cycle <= 0.34) {
    u = 35 - (35 / 0.02) * (t_in_cycle - 0.32);
  } else if (t_in_cycle > 0.34 && t_in_cycle <= 0.46) {
    u = 24 - 12 / 0.12 * (t_in_cycle-0.1);
  } else if (t_in_cycle > 0.46 && t_in_cycle <= 0.59) {
    u = -12;
  } else if (t_in_cycle > 0.59 && t_in_cycle <= 0.6) {
    u = -443.2 + 8.8 / 0.01 * (t_in_cycle-0.1);
  } else {
    u = -3.2;
  }

  if (u >= 0.0) {
    u_plus = u;
    u_minus = 0.0;
  } else {
    u_plus = 0.0;
    u_minus = u;
  }
}
