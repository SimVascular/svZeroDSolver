// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ClosedLoopHeartPulmonary.h"

#include "Model.h"

void ClosedLoopHeartPulmonary::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 14,
                     {"V_RA", "Q_RA", "P_RV", "V_RV", "Q_RV", "P_pul", "P_LA",
                      "V_LA", "Q_LA", "P_LV", "V_LV", "Q_LV"});
}

void ClosedLoopHeartPulmonary::update_constant(
    SparseSystem &system, std::vector<double> &parameters) {
  // DOF 0, Eq 0: Right atrium pressure
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;

  // DOF 2, Eq 1: Aortic pressure
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[2]) =
      parameters[global_param_ids[ParamId::CPA]];
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = 1.0;

  // DOF 4, Eq 2: Right atrium volume
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[4]) = 1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[1]) = -1.0;

  // DOF 5, Eq 3: Right atrium outflow
  system.E.coeffRef(global_eqn_ids[3], global_var_ids[5]) =
      parameters[global_param_ids[ParamId::LRA_V]];
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[0]) = -1.0;
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[6]) = 1.0;

  // DOF 6, Eq 4: Right ventricle pressure
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[6]) = 1.0;

  // DOF 7, Eq 5: Right ventricle volume
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[7]) = 1.0;

  // DOF 8, Eq 6: Right ventricle outflow
  system.E.coeffRef(global_eqn_ids[6], global_var_ids[8]) =
      parameters[global_param_ids[ParamId::LRV_A]];
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[6]) = -1.0;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[9]) = 1.0;

  // DOF 9, Eq 7: Pulmonary pressure
  system.E.coeffRef(global_eqn_ids[7], global_var_ids[9]) =
      parameters[global_param_ids[ParamId::CP]];
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[9]) =
      1.0 / parameters[global_param_ids[ParamId::RPD]];
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[10]) =
      -1.0 / parameters[global_param_ids[ParamId::RPD]];

  // DOF 10, Eq 8: Left atrium pressure
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[10]) = 1.0;

  // DOF 11, Eq 9: Left atrium volume
  system.E.coeffRef(global_eqn_ids[9], global_var_ids[11]) = 1.0;

  // DOF 12, Eq 10: Left atrium outflow
  system.E.coeffRef(global_eqn_ids[10], global_var_ids[12]) =
      parameters[global_param_ids[ParamId::LLA_V]];
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[10]) = -1.0;
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[13]) = 1.0;

  // DOF 13, Eq 11: Left ventricle pressure
  system.F.coeffRef(global_eqn_ids[11], global_var_ids[13]) = 1.0;

  // DOF 14, Eq 12: Left ventricle volume
  system.E.coeffRef(global_eqn_ids[12], global_var_ids[14]) = 1.0;

  // DOF 15, Eq 13: Left ventricle outflow
  system.F.coeffRef(global_eqn_ids[13], global_var_ids[2]) = 1.0;
  system.F.coeffRef(global_eqn_ids[13], global_var_ids[13]) = -1.0;
  system.E.coeffRef(global_eqn_ids[13], global_var_ids[15]) =
      parameters[global_param_ids[ParamId::LLV_A]];
}

void ClosedLoopHeartPulmonary::update_time(SparseSystem &system,
                                           std::vector<double> &parameters) {
  get_activation_and_elastance_functions(parameters);

  // DOF 0, Eq 0: Right atrium pressure
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      -AA * parameters[global_param_ids[ParamId::EMAX_RA]];

  // DOF 6, Eq 4: Right ventricle pressure
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[7]) = -Erv;
  system.C(global_eqn_ids[4]) =
      Erv * parameters[global_param_ids[ParamId::VRV_U]];

  // DOF 10, Eq 8: Left atrium pressure
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[11]) =
      -AA * parameters[global_param_ids[ParamId::EMAX_LA]];

  // DOF 13, Eq 11: Left ventricle pressure
  system.F.coeffRef(global_eqn_ids[11], global_var_ids[14]) = -Elv;
  system.C(global_eqn_ids[11]) =
      Elv * parameters[global_param_ids[ParamId::VLV_U]];
}

void ClosedLoopHeartPulmonary::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {
  get_psi_ra_la(parameters, y);
  get_valve_positions(y);

  // Technically, F matrix and C vector neither depend on time nor solution
  // However, we treat F here as constant (despite the solution-dependent
  // update) since only 0/1 entries of the valves change, which are assumed
  // constant in the linearization. Thus, F behaves here like a constant block
  // for the assembly in SparseSystem

  // DOF IDs are arranged as inflow
  // [P_in,Q_in,P_out,Q_out,internal variables...]

  // DOF 0, Eq 0: Right atrium pressure
  system.C(global_eqn_ids[0]) =
      AA * parameters[global_param_ids[ParamId::EMAX_RA]] *
          parameters[global_param_ids[ParamId::VASO_RA]] +
      psi_ra * (AA - 1.0);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      psi_ra_derivative * (AA - 1.0);

  // DOF 10, Eq 8: Left atrium pressure
  system.C(global_eqn_ids[8]) =
      AA * parameters[global_param_ids[ParamId::EMAX_LA]] *
          parameters[global_param_ids[ParamId::VASO_LA]] +
      psi_la * (AA - 1.0);
  system.dC_dy.coeffRef(global_eqn_ids[8], global_var_ids[11]) =
      psi_la_derivative * (AA - 1.0);

  // DOF 2, Eq 1: Aortic pressure
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[15]) = -valves[15];

  // DOF 9, Eq 7: Pulmonary pressure
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[8]) = -valves[8];

  // DOF 4, Eq 2: Right atrium volume
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[5]) = valves[5];

  // DOF 7, Eq 5: Right ventricle volume
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[5]) = -valves[5];
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[8]) = valves[8];

  // DOF 11, Eq 9: Left atrium volume
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[8]) = -valves[8];
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[12]) = valves[12];

  // DOF 14, Eq 12: Left ventricle volume
  system.F.coeffRef(global_eqn_ids[12], global_var_ids[12]) = -valves[12];
  system.F.coeffRef(global_eqn_ids[12], global_var_ids[15]) = valves[15];

  // DOF 5, Eq 3: Right atrium outflow
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[5]) =
      parameters[global_param_ids[ParamId::RRA_V]] * valves[5];

  // DOF 8, Eq 6: Right ventricle outflow
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[8]) =
      parameters[global_param_ids[ParamId::RRV_A]] * valves[8];

  // DOF 12, Eq 10: Left atrium outflow
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[12]) =
      parameters[global_param_ids[ParamId::RLA_V]] * valves[12];

  // DOF 15, Eq 13: Left ventricle outflow
  system.F.coeffRef(global_eqn_ids[13], global_var_ids[15]) =
      parameters[global_param_ids[ParamId::RLV_AO]] * valves[15];
}

void ClosedLoopHeartPulmonary::get_activation_and_elastance_functions(
    std::vector<double> &parameters) {
  auto T_cardiac = model->cardiac_cycle_period;
  auto Tsa = T_cardiac * parameters[global_param_ids[ParamId::TSA]];
  auto tpwave = T_cardiac / parameters[global_param_ids[ParamId::TPWAVE]];
  auto t_in_cycle = fmod(model->time, T_cardiac);

  // Activation function
  AA = 0.0;
  if (t_in_cycle <= tpwave) {
    AA = (0.5) * (1.0 - cos(2.0 * M_PI * (t_in_cycle - tpwave + Tsa) / Tsa));
  } else if ((t_in_cycle >= (T_cardiac - Tsa) + tpwave) &&
             (t_in_cycle < T_cardiac)) {
    AA = (0.5) * (1.0 - cos(2.0 * M_PI *
                            (t_in_cycle - tpwave - (T_cardiac - Tsa)) / Tsa));
  } else {
    AA = 0.0;
  }

  // Elastance modes (copied from J. Tran's tuning framework)
  const int num_elast_modes = 25;
  double Ft_elastance[num_elast_modes][2] = {
      {0.283748803, 0.000000000},   {0.031830626, -0.374299825},
      {-0.209472400, -0.018127770}, {0.020520047, 0.073971113},
      {0.008316883, -0.047249597},  {-0.041677660, 0.003212163},
      {0.000867323, 0.019441411},   {-0.001675379, -0.005565534},
      {-0.011252277, 0.003401432},  {-0.000414677, 0.008376795},
      {0.000253749, -0.000071880},  {-0.002584966, 0.001566861},
      {0.000584752, 0.003143555},   {0.000028502, -0.000024787},
      {0.000022961, -0.000007476},  {0.000018735, -0.000001281},
      {0.000015573, 0.000001781},   {0.000013133, 0.000003494},
      {0.000011199, 0.000004507},   {0.000009634, 0.000005117},
      {0.000008343, 0.000005481},   {0.000007265, 0.000005687},
      {0.000006354, 0.000005789},   {0.000005575, 0.000005821},
      {0.000004903, 0.000005805}};

  // RV and LV elastance
  double Elv_i = 0.0;

  for (auto i = 0; i < num_elast_modes; i++) {
    Elv_i =
        Elv_i +
        (Ft_elastance[i][0]) * cos(2.0 * M_PI * i * t_in_cycle / T_cardiac) -
        (Ft_elastance[i][1]) * sin(2.0 * M_PI * i * t_in_cycle / T_cardiac);
  }

  Elv = Elv_i * parameters[global_param_ids[ParamId::ELV_S]];
  Erv = Elv_i * parameters[global_param_ids[ParamId::ERV_S]];
}

void ClosedLoopHeartPulmonary::get_psi_ra_la(
    std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  auto RA_volume = y[global_var_ids[4]];
  auto LA_volume = y[global_var_ids[11]];
  psi_ra = parameters[global_param_ids[ParamId::KXP_RA]] *
           (exp((RA_volume - parameters[global_param_ids[ParamId::VASO_RA]]) *
                parameters[global_param_ids[ParamId::KXV_RA]]) -
            1.0);
  psi_la = parameters[global_param_ids[ParamId::KXP_LA]] *
           (exp((LA_volume - parameters[global_param_ids[ParamId::VASO_LA]]) *
                parameters[global_param_ids[ParamId::KXV_LA]]) -
            1.0);

  psi_ra_derivative =
      parameters[global_param_ids[ParamId::KXP_RA]] *
      exp((RA_volume - parameters[global_param_ids[ParamId::VASO_RA]]) *
          parameters[global_param_ids[ParamId::KXV_RA]]) *
      parameters[global_param_ids[ParamId::KXV_RA]];
  psi_la_derivative =
      parameters[global_param_ids[ParamId::KXP_LA]] *
      exp((LA_volume - parameters[global_param_ids[ParamId::VASO_LA]]) *
          parameters[global_param_ids[ParamId::KXV_LA]]) *
      parameters[global_param_ids[ParamId::KXV_LA]];
}

void ClosedLoopHeartPulmonary::get_valve_positions(
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  std::fill(valves, valves + 16, 1.0);

  // RA to RV
  auto pressure_ra = y[global_var_ids[0]];
  auto pressure_rv = y[global_var_ids[6]];
  auto outflow_ra = y[global_var_ids[5]];
  if ((pressure_ra <= pressure_rv) && (outflow_ra <= 0.0)) {
    valves[5] = 0.0;
  }

  // RV to pulmonary
  auto pressure_pulmonary = y[global_var_ids[9]];
  auto outflow_rv = y[global_var_ids[8]];
  if ((pressure_rv <= pressure_pulmonary) && (outflow_rv <= 0.0)) {
    valves[8] = 0.0;
  }

  // LA to LV
  auto pressure_la = y[global_var_ids[10]];
  auto pressure_lv = y[global_var_ids[13]];
  auto outflow_la = y[global_var_ids[12]];
  if ((pressure_la <= pressure_lv) && (outflow_la <= 0.0)) {
    valves[12] = 0.0;
  }

  // LV to aorta
  auto pressure_aorta = y[global_var_ids[2]];
  auto outflow_lv = y[global_var_ids[15]];
  if ((pressure_lv <= pressure_aorta) && (outflow_lv <= 0.0)) {
    valves[15] = 0.0;
  }
}

void ClosedLoopHeartPulmonary::post_solve(
    Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  for (size_t i = 0; i < 16; i++)
    if (valves[i] < 0.5) y[global_var_ids[i]] = 0.0;
}
