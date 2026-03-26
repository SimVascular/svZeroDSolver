// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ClosedLoopHeartPulmonarySmooth.h"
#include "Model.h"

// Identical to ClosedLoopHeartPulmonary except get_valve_positions uses
// smooth tanh and there is no post_solve.

void ClosedLoopHeartPulmonarySmooth::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 14,
                     {"V_RA", "Q_RA", "P_RV", "V_RV", "Q_RV", "P_pul", "P_LA",
                      "V_LA", "Q_LA", "P_LV", "V_LV", "Q_LV"});
}

// update_constant: IDENTICAL to original
void ClosedLoopHeartPulmonarySmooth::update_constant(
    SparseSystem& system, std::vector<double>& parameters) {
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[2]) = parameters[global_param_ids[CPA]];
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = 1.0;
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[4]) = 1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[1]) = -1.0;
  system.E.coeffRef(global_eqn_ids[3], global_var_ids[5]) = parameters[global_param_ids[LRA_V]];
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[0]) = -1.0;
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[6]) = 1.0;
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[6]) = 1.0;
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[7]) = 1.0;
  system.E.coeffRef(global_eqn_ids[6], global_var_ids[8]) = parameters[global_param_ids[LRV_A]];
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[6]) = -1.0;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[9]) = 1.0;
  system.E.coeffRef(global_eqn_ids[7], global_var_ids[9]) = parameters[global_param_ids[CP]];
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[9]) = 1.0 / parameters[global_param_ids[RPD]];
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[10]) = -1.0 / parameters[global_param_ids[RPD]];
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[10]) = 1.0;
  system.E.coeffRef(global_eqn_ids[9], global_var_ids[11]) = 1.0;
  system.E.coeffRef(global_eqn_ids[10], global_var_ids[12]) = parameters[global_param_ids[LLA_V]];
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[10]) = -1.0;
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[13]) = 1.0;
  system.F.coeffRef(global_eqn_ids[11], global_var_ids[13]) = 1.0;
  system.E.coeffRef(global_eqn_ids[12], global_var_ids[14]) = 1.0;
  system.F.coeffRef(global_eqn_ids[13], global_var_ids[2]) = 1.0;
  system.F.coeffRef(global_eqn_ids[13], global_var_ids[13]) = -1.0;
  system.E.coeffRef(global_eqn_ids[13], global_var_ids[15]) = parameters[global_param_ids[LLV_A]];
}

// update_time: IDENTICAL to original
void ClosedLoopHeartPulmonarySmooth::update_time(
    SparseSystem& system, std::vector<double>& parameters) {
  get_activation_and_elastance_functions(parameters);
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -AA * parameters[global_param_ids[EMAX_RA]];
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[7]) = -Erv;
  system.C(global_eqn_ids[4]) = Erv * parameters[global_param_ids[VRV_U]];
  system.F.coeffRef(global_eqn_ids[8], global_var_ids[11]) = -AA * parameters[global_param_ids[EMAX_LA]];
  system.F.coeffRef(global_eqn_ids[11], global_var_ids[14]) = -Elv;
  system.C(global_eqn_ids[11]) = Elv * parameters[global_param_ids[VLV_U]];
}

// update_solution: IDENTICAL to original (valve-gated F entries treated as constant)
void ClosedLoopHeartPulmonarySmooth::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  get_psi_ra_la(parameters, y);

  // Evaluate valve state once per time step (frozen during Newton iterations)
  // for correct Jacobian (F treated as constant)
  if (model->time != last_time_) {
    last_time_ = model->time;
    get_valve_positions(parameters, y);
  }

  // Atrial P-V
  system.C(global_eqn_ids[0]) = AA * parameters[global_param_ids[EMAX_RA]] * parameters[global_param_ids[VASO_RA]] + psi_ra * (AA - 1.0);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) = psi_ra_derivative * (AA - 1.0);
  system.C(global_eqn_ids[8]) = AA * parameters[global_param_ids[EMAX_LA]] * parameters[global_param_ids[VASO_LA]] + psi_la * (AA - 1.0);
  system.dC_dy.coeffRef(global_eqn_ids[8], global_var_ids[11]) = psi_la_derivative * (AA - 1.0);

  // Valve sigmoid derivatives: ds/dP_up = k*s*(1-s), ds/dP_down = -k*s*(1-s)
  double k = parameters[global_param_ids[STEEPNESS]];
  double s5=valves[5], s8=valves[8], s12=valves[12], s15=valves[15];
  double d5=k*s5*(1-s5), d8=k*s8*(1-s8), d12=k*s12*(1-s12), d15=k*s15*(1-s15);
  double Qa=y[global_var_ids[5]], Qr=y[global_var_ids[8]], Ql=y[global_var_ids[12]], Qv=y[global_var_ids[15]];
  double Ra=parameters[global_param_ids[RRA_V]], Rr=parameters[global_param_ids[RRV_A]];
  double Rl=parameters[global_param_ids[RLA_V]], Rv=parameters[global_param_ids[RLV_AO]];

  // Valve-gated terms in F (treated as constant during linearization)
  double Rmax = parameters[global_param_ids[RMAX]];
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[15]) = -s15;
  system.F.coeffRef(global_eqn_ids[7], global_var_ids[8]) = -s8;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[5]) = s5;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[5]) = -s5;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[8]) = s8;
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[8]) = -s8;
  system.F.coeffRef(global_eqn_ids[9], global_var_ids[12]) = s12;
  system.F.coeffRef(global_eqn_ids[12], global_var_ids[12]) = -s12;
  system.F.coeffRef(global_eqn_ids[12], global_var_ids[15]) = s15;
  // Outflow: R_eff = R*s + Rmax*(1-s)
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[5]) = Ra*s5+Rmax*(1-s5);
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[8]) = Rr*s8+Rmax*(1-s8);
  system.F.coeffRef(global_eqn_ids[10], global_var_ids[12]) = Rl*s12+Rmax*(1-s12);
  system.F.coeffRef(global_eqn_ids[13], global_var_ids[15]) = Rv*s15+Rmax*(1-s15);
}

// Activation/elastance: IDENTICAL to original
void ClosedLoopHeartPulmonarySmooth::get_activation_and_elastance_functions(std::vector<double>& parameters) {
  auto T = model->cardiac_cycle_period;
  auto Tsa = T * parameters[global_param_ids[TSA]];
  auto tpw = T / parameters[global_param_ids[TPWAVE]];
  auto tc = fmod(model->time, T);
  AA = 0.0;
  if (tc <= tpw) AA = 0.5*(1.0-cos(2.0*M_PI*(tc-tpw+Tsa)/Tsa));
  else if (tc >= (T-Tsa)+tpw && tc < T) AA = 0.5*(1.0-cos(2.0*M_PI*(tc-tpw-(T-Tsa))/Tsa));
  const int N=25;
  double Ft[N][2]={{0.283748803,0},{.031830626,-.374299825},{-.2094724,-.01812777},{.020520047,.073971113},{.008316883,-.047249597},{-.04167766,.003212163},{.000867323,.019441411},{-.001675379,-.005565534},{-.011252277,.003401432},{-.000414677,.008376795},{.000253749,-.00007188},{-.002584966,.001566861},{.000584752,.003143555},{.000028502,-.000024787},{.000022961,-.000007476},{.000018735,-.000001281},{.000015573,.000001781},{.000013133,.000003494},{.000011199,.000004507},{.000009634,.000005117},{.000008343,.000005481},{.000007265,.000005687},{.000006354,.000005789},{.000005575,.000005821},{.000004903,.000005805}};
  double Ei=0;
  for(int i=0;i<N;i++) Ei+=Ft[i][0]*cos(2*M_PI*i*tc/T)-Ft[i][1]*sin(2*M_PI*i*tc/T);
  Elv=Ei*parameters[global_param_ids[ELV_S]];
  Erv=Ei*parameters[global_param_ids[ERV_S]];
}

void ClosedLoopHeartPulmonarySmooth::get_psi_ra_la(std::vector<double>& parameters, const Eigen::Matrix<double,Eigen::Dynamic,1>& y) {
  auto Vra=y[global_var_ids[4]], Vla=y[global_var_ids[11]];
  double kp_r=parameters[global_param_ids[KXP_RA]],kv_r=parameters[global_param_ids[KXV_RA]],v_r=parameters[global_param_ids[VASO_RA]];
  double kp_l=parameters[global_param_ids[KXP_LA]],kv_l=parameters[global_param_ids[KXV_LA]],v_l=parameters[global_param_ids[VASO_LA]];
  psi_ra=kp_r*(exp((Vra-v_r)*kv_r)-1); psi_la=kp_l*(exp((Vla-v_l)*kv_l)-1);
  psi_ra_derivative=kp_r*exp((Vra-v_r)*kv_r)*kv_r;
  psi_la_derivative=kp_l*exp((Vla-v_l)*kv_l)*kv_l;
}

// THE ONLY DIFFERENCE: smooth tanh valves, no post_solve
void ClosedLoopHeartPulmonarySmooth::get_valve_positions(std::vector<double>& parameters, const Eigen::Matrix<double,Eigen::Dynamic,1>& y) {
  double k = parameters[global_param_ids[STEEPNESS]];
  valves[5]  = 0.5*(1.0+tanh(k*(y[global_var_ids[0]]-y[global_var_ids[6]])));
  valves[8]  = 0.5*(1.0+tanh(k*(y[global_var_ids[6]]-y[global_var_ids[9]])));
  valves[12] = 0.5*(1.0+tanh(k*(y[global_var_ids[10]]-y[global_var_ids[13]])));
  valves[15] = 0.5*(1.0+tanh(k*(y[global_var_ids[13]]-y[global_var_ids[2]])));
}
