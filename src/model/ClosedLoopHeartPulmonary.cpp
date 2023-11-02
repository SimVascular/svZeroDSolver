// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "ClosedLoopHeartPulmonary.h"

#include "Model.h"

void ClosedLoopHeartPulmonary::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 14,
                     {"V_RA", "Q_RA", "P_RV", "V_RV", "Q_RV", "P_pul", "P_LA",
                      "V_LA", "Q_LA", "P_LV", "V_LV", "Q_LV"});
}

void ClosedLoopHeartPulmonary::update_constant(
    SparseSystem &system, std::vector<double> &parameters) {
  // DOF 2, Eq 1: Aortic pressure
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) =
      parameters[this->global_param_ids[ParamId::CPA]];
  // DOF 4, Eq 2: Right atrium volume
  system.E.coeffRef(this->global_eqn_ids[2], this->global_var_ids[4]) = 1.0;
  // DOF 5, Eq 3: Right atrium outflow
  system.E.coeffRef(this->global_eqn_ids[3], this->global_var_ids[5]) =
      parameters[this->global_param_ids[ParamId::LRA_V]];
  // DOF 7, Eq 5: Right ventricle volume
  system.E.coeffRef(this->global_eqn_ids[5], this->global_var_ids[7]) = 1.0;
  // DOF 8, Eq 6: Right ventricle outflow
  system.E.coeffRef(this->global_eqn_ids[6], this->global_var_ids[8]) =
      parameters[this->global_param_ids[ParamId::LRV_A]];
  // DOF 9, Eq 7: Pulmonary pressure
  system.E.coeffRef(this->global_eqn_ids[7], this->global_var_ids[9]) =
      parameters[this->global_param_ids[ParamId::CP]];
  // DOF 11, Eq 9: Left atrium volume
  system.E.coeffRef(this->global_eqn_ids[9], this->global_var_ids[11]) = 1.0;
  // DOF 12, Eq 10: Left atrium outflow
  system.E.coeffRef(this->global_eqn_ids[10], this->global_var_ids[12]) =
      parameters[this->global_param_ids[ParamId::LLA_V]];
  // DOF 14, Eq 12: Left ventricle volume
  system.E.coeffRef(this->global_eqn_ids[12], this->global_var_ids[14]) = 1.0;
  // DOF 15, Eq 13: Left ventricle outflow
  system.E.coeffRef(this->global_eqn_ids[13], this->global_var_ids[15]) =
      parameters[this->global_param_ids[ParamId::LLV_A]];
}

void ClosedLoopHeartPulmonary::update_time(SparseSystem &system,
                                           std::vector<double> &parameters) {
  this->get_activation_and_elastance_functions(parameters);
}

void ClosedLoopHeartPulmonary::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {
  this->get_psi_ra_la(parameters, y);
  this->get_valve_positions(y);

  // F and C matrices depend on time and solution
  // Specifying all terms here, including constant terms (which can instead be
  // specified in update_constant) for readability (Doesn't seem to make a
  // difference to compute time) DOF IDs are arranged as inflow
  // [P_in,Q_in,P_out,Q_out,internal variables...]

  // DOF 0, Eq 0: Right atrium pressure
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[4]) =
      -this->AA * parameters[this->global_param_ids[ParamId::EMAX_RA]];
  system.C(this->global_eqn_ids[0]) =
      this->AA * parameters[this->global_param_ids[ParamId::EMAX_RA]] *
          parameters[this->global_param_ids[ParamId::VASO_RA]] +
      psi_ra * (this->AA - 1.0);
  system.D.coeffRef(this->global_eqn_ids[0], this->global_var_ids[4]) =
      psi_ra_derivative * (this->AA - 1.0);

  // DOF 1: Flow into right atrium (no equation)

  // DOF 2, Eq 1: Aortic pressure
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[15]) =
      -valves[15];
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[3]) = 1.0;

  // DOF 3: Flow into aorta (no equation)

  // DOF 4, Eq 2: Right atrium volume
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[5]) =
      1.0 * valves[5];
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[1]) = -1.0;

  // DOF 5, Eq 3: Right atrium outflow
  system.F.coeffRef(this->global_eqn_ids[3], this->global_var_ids[5]) =
      parameters[this->global_param_ids[ParamId::RRA_V]] * valves[5];
  system.F.coeffRef(this->global_eqn_ids[3], this->global_var_ids[0]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[3], this->global_var_ids[6]) = 1.0;

  // DOF 6, Eq 4: Right ventricle pressure
  system.F.coeffRef(this->global_eqn_ids[4], this->global_var_ids[6]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[4], this->global_var_ids[7]) =
      -this->Erv;
  system.C(this->global_eqn_ids[4]) =
      this->Erv * parameters[this->global_param_ids[ParamId::VRV_U]];

  // DOF 7, Eq 5: Right ventricle volume
  system.F.coeffRef(this->global_eqn_ids[5], this->global_var_ids[5]) =
      -1.0 * valves[5];
  system.F.coeffRef(this->global_eqn_ids[5], this->global_var_ids[8]) =
      1.0 * valves[8];

  // DOF 8, Eq 6: Right ventricle outflow
  system.F.coeffRef(this->global_eqn_ids[6], this->global_var_ids[6]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[6], this->global_var_ids[9]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[6], this->global_var_ids[8]) =
      parameters[this->global_param_ids[ParamId::RRV_A]] * valves[8];

  // DOF 9, Eq 7: Pulmonary pressure
  system.F.coeffRef(this->global_eqn_ids[7], this->global_var_ids[8]) =
      -valves[8];
  system.F.coeffRef(this->global_eqn_ids[7], this->global_var_ids[9]) =
      1.0 / parameters[this->global_param_ids[ParamId::RPD]];
  system.F.coeffRef(this->global_eqn_ids[7], this->global_var_ids[10]) =
      -1.0 / parameters[this->global_param_ids[ParamId::RPD]];

  // DOF 10, Eq 8: Left atrium pressure
  system.F.coeffRef(this->global_eqn_ids[8], this->global_var_ids[10]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[8], this->global_var_ids[11]) =
      -this->AA * parameters[this->global_param_ids[ParamId::EMAX_LA]];
  system.C(this->global_eqn_ids[8]) =
      this->AA * parameters[this->global_param_ids[ParamId::EMAX_LA]] *
          parameters[this->global_param_ids[ParamId::VASO_LA]] +
      psi_la * (this->AA - 1.0);
  system.D.coeffRef(this->global_eqn_ids[8], this->global_var_ids[11]) =
      psi_la_derivative * (this->AA - 1.0);

  // DOF 11, Eq 9: Left atrium volume
  system.F.coeffRef(this->global_eqn_ids[9], this->global_var_ids[8]) =
      -1.0 * valves[8];
  system.F.coeffRef(this->global_eqn_ids[9], this->global_var_ids[12]) =
      1.0 * valves[12];

  // DOF 12, Eq 10: Left atrium outflow
  system.F.coeffRef(this->global_eqn_ids[10], this->global_var_ids[10]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[10], this->global_var_ids[13]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[10], this->global_var_ids[12]) =
      parameters[this->global_param_ids[ParamId::RLA_V]] * valves[12];

  // DOF 13, Eq 11: Left ventricle pressure
  system.F.coeffRef(this->global_eqn_ids[11], this->global_var_ids[13]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[11], this->global_var_ids[14]) =
      -this->Elv;
  system.C(this->global_eqn_ids[11]) =
      this->Elv * parameters[this->global_param_ids[ParamId::VLV_U]];

  // DOF 14, Eq 12: Left ventricle volume
  system.F.coeffRef(this->global_eqn_ids[12], this->global_var_ids[12]) =
      -1.0 * valves[12];
  system.F.coeffRef(this->global_eqn_ids[12], this->global_var_ids[15]) =
      1.0 * valves[15];

  // DOF 15, Eq 13: Left ventricle outflow
  system.F.coeffRef(this->global_eqn_ids[13], this->global_var_ids[13]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[13], this->global_var_ids[2]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[13], this->global_var_ids[15]) =
      parameters[this->global_param_ids[ParamId::RLV_AO]] * valves[15];
}

void ClosedLoopHeartPulmonary::get_activation_and_elastance_functions(
    std::vector<double> &parameters) {
  auto T_cardiac = this->model->cardiac_cycle_period;
  auto Tsa = T_cardiac * parameters[this->global_param_ids[ParamId::TSA]];
  auto tpwave = T_cardiac / parameters[this->global_param_ids[ParamId::TPWAVE]];
  auto t_in_cycle = fmod(this->model->time, T_cardiac);

  // Activation function
  AA = 0.0;
  if (t_in_cycle <= tpwave) {
    AA = (0.5) * (1.0 - cos(2.0 * PI * (t_in_cycle - tpwave + Tsa) / Tsa));
  } else if ((t_in_cycle >= (T_cardiac - Tsa) + tpwave) and
             (t_in_cycle < T_cardiac)) {
    AA =
        (0.5) *
        (1.0 - cos(2.0 * PI * (t_in_cycle - tpwave - (T_cardiac - Tsa)) / Tsa));
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
    Elv_i = Elv_i +
            (Ft_elastance[i][0]) * cos(2.0 * PI * i * t_in_cycle / T_cardiac) -
            (Ft_elastance[i][1]) * sin(2.0 * PI * i * t_in_cycle / T_cardiac);
  }

  Elv = Elv_i * parameters[this->global_param_ids[ParamId::ELV_S]];
  Erv = Elv_i * parameters[this->global_param_ids[ParamId::ERV_S]];
}

void ClosedLoopHeartPulmonary::get_psi_ra_la(
    std::vector<double> &parameters,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  auto RA_volume = y[this->global_var_ids[4]];
  auto LA_volume = y[this->global_var_ids[11]];
  psi_ra =
      parameters[this->global_param_ids[ParamId::KXP_RA]] *
      (exp((RA_volume - parameters[this->global_param_ids[ParamId::VASO_RA]]) *
           parameters[this->global_param_ids[ParamId::KXV_RA]]) -
       1.0);
  psi_la =
      parameters[this->global_param_ids[ParamId::KXP_LA]] *
      (exp((LA_volume - parameters[this->global_param_ids[ParamId::VASO_LA]]) *
           parameters[this->global_param_ids[ParamId::KXV_LA]]) -
       1.0);

  psi_ra_derivative =
      parameters[this->global_param_ids[ParamId::KXP_RA]] *
      exp((RA_volume - parameters[this->global_param_ids[ParamId::VASO_RA]]) *
          parameters[this->global_param_ids[ParamId::KXV_RA]]) *
      parameters[this->global_param_ids[ParamId::KXV_RA]];
  psi_la_derivative =
      parameters[this->global_param_ids[ParamId::KXP_LA]] *
      exp((LA_volume - parameters[this->global_param_ids[ParamId::VASO_LA]]) *
          parameters[this->global_param_ids[ParamId::KXV_LA]]) *
      parameters[this->global_param_ids[ParamId::KXV_LA]];
}

void ClosedLoopHeartPulmonary::get_valve_positions(
    Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  std::fill(valves, valves + 16, 1.0);

  // RA to RV
  auto pressure_ra = y[this->global_var_ids[0]];
  auto pressure_rv = y[this->global_var_ids[6]];
  auto outflow_ra = y[this->global_var_ids[5]];
  if ((pressure_ra <= pressure_rv) and (outflow_ra <= 0.0)) {
    valves[5] = 0.0;
    y[this->global_var_ids[5]] = 0.0;
  }

  // RV to pulmonary
  auto pressure_pulmonary = y[this->global_var_ids[9]];
  auto outflow_rv = y[this->global_var_ids[8]];
  if ((pressure_rv <= pressure_pulmonary) and (outflow_rv <= 0.0)) {
    valves[8] = 0.0;
    y[this->global_var_ids[8]] = 0.0;
  }

  // LA to LV
  auto pressure_la = y[this->global_var_ids[10]];
  auto pressure_lv = y[this->global_var_ids[13]];
  auto outflow_la = y[this->global_var_ids[12]];
  if ((pressure_la <= pressure_lv) and (outflow_la <= 0.0)) {
    valves[12] = 0.0;
    y[this->global_var_ids[12]] = 0.0;
  }

  // LV to aorta
  auto pressure_aorta = y[this->global_var_ids[2]];
  auto outflow_lv = y[this->global_var_ids[15]];
  if ((pressure_lv <= pressure_aorta) and (outflow_lv <= 0.0)) {
    valves[15] = 0.0;
    y[this->global_var_ids[15]] = 0.0;
  }
}

