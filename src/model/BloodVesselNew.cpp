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

#include "BloodVesselNew.h"

void BloodVesselNew::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 8, {"r","v","S","k","tau","V"});
}

void BloodVesselNew::update_constant(SparseSystem &system,
                                  std::vector<double> &parameters) {
  // Hard code constant parameters

  // const double rho = 1000.0; 
  // const double d = 0.0142; 
  // const double Ro = 0.0236; 
  // const double W1 = 1.0000; 
  // const double W2 = 1.0000; 
  // const double eta = 70.0000; 
  // const double a = 1.0000; 
  // const double sigma_o = 124000.0000; 


   // Set element contributions
  // coeffRef args are the indices (i,j) of the matrix
  // global_eqn_ids: number of rows in the matrix, set in setup_dofs
  // global_var_ids: number of columns, organized as pressure and flow of all
  // inlets and then all outlets of the block
 
// deleted E and F matrices %%%%%%
}

// void ChamberElastanceInductor::update_time(SparseSystem &system,
//   std::vector<double> &parameters) {
// get_elastance_values(parameters);

// system.C.coeffRef(global_eqn_ids[2]) = tau_ - a * sigma_o + a * tau;
// system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[8]) = a;
// }

// void ChamberElastanceInductor::get_elastance_values(
// std::vector<double> &parameters) {
// double Emax = parameters[global_param_ids[ParamId::EMAX]];
// double Emin = parameters[global_param_ids[ParamId::EMIN]];
// double Vrd = parameters[global_param_ids[ParamId::VRD]];
// double Vrs = parameters[global_param_ids[ParamId::VRS]];
// double t_active = parameters[global_param_ids[ParamId::TACTIVE]];
// double t_twitch = parameters[global_param_ids[ParamId::TTWITCH]];

// auto T_cardiac = model->cardiac_cycle_period;
// auto t_in_cycle = fmod(model->time, T_cardiac);

// double t_contract = 0;
// if (t_in_cycle >= t_active) {
// t_contract = t_in_cycle - t_active;
// }

// double act = 0;
// if (t_contract <= t_twitch) {
// act = -0.5 * cos(2 * M_PI * t_contract / t_twitch) + 0.5;
// }

// Vrest = (1.0 - act) * (Vrd - Vrs) + Vrs;
// a = (Emax - Emin) * act + Emin; // CHANGED "ELAS" TO "a" TO MATCH REST OF SCRIPT

// }


void BloodVesselNew::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {
  // Get parameters
  double rho = parameters[global_param_ids[ParamId::rho]];
  std::cout<<rho<<std::endl;
  double d = parameters[global_param_ids[ParamId::d]];
  double Ro = parameters[global_param_ids[ParamId::Ro]];
  double W1 = parameters[global_param_ids[ParamId::W1]];
  double W2 = parameters[global_param_ids[ParamId::W2]];
  double eta = parameters[global_param_ids[ParamId::eta]];
  double a = parameters[global_param_ids[ParamId::a]];
  double sigma_o = parameters[global_param_ids[ParamId::sigma_o]];

  double Pin = y[global_var_ids[0]]; 
  double Qin = y[global_var_ids[1]]; 
  double Pout = y[global_var_ids[2]]; 
  double Qout = y[global_var_ids[3]]; 
  double r = y[global_var_ids[4]]; 
  double v = y[global_var_ids[5]]; 
  double S = y[global_var_ids[6]]; 
  double k = y[global_var_ids[7]]; 
  double tau = y[global_var_ids[8]]; 
  double V = y[global_var_ids[9]];

  double Pin_ = dy[global_var_ids[0]]; 
  double Qin_ = dy[global_var_ids[1]]; 
  double Pout_ = dy[global_var_ids[2]]; 
  double Qout_ = dy[global_var_ids[3]]; 
  double r_ = dy[global_var_ids[4]]; 
  double v_ = dy[global_var_ids[5]]; 
  double S_ = dy[global_var_ids[6]]; 
  double k_ = dy[global_var_ids[7]]; 
  double tau_ = dy[global_var_ids[8]]; 
  double V_ = dy[global_var_ids[9]];
  // act_stress = get_active_stress(parameters); %%%%%% come back to add active component
  

  // Set element contributions
  system.C(global_eqn_ids[0]) = -Pout*pow(r/Ro+1.0,2.0)+d*rho*v_+(S*d*(r/Ro+1.0))/Ro;
  system.C(global_eqn_ids[1]) = -S+tau-(1.0/(k*k*k)*4.0-4.0)*(W1+W2*k)-eta*k_*(1.0/(k*k*k*k*k*k)*2.0-1.0)*2.0;
  system.C(global_eqn_ids[2]) = tau_-a*sigma_o+a*tau;
  system.C(global_eqn_ids[3]) = Qin-Qout-V_;
  system.C(global_eqn_ids[4]) = -V_+(Ro*Ro)*v*M_PI*pow(r/Ro+1.0,2.0)*4.0;
  system.C(global_eqn_ids[5]) = k-pow(r/Ro+1.0,2.0);
  system.C(global_eqn_ids[6]) = r_-v;
  system.C(global_eqn_ids[7]) = Pin-Pout;

  std::cout<<"wtf1"<<std::endl;

  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -pow(r/Ro+1.0,2.0);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -1.0/(Ro*Ro)*(Pout*Ro*2.0-S*d+Pout*r*2.0);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[6]) = 1.0/(Ro*Ro)*d*(Ro+r);
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[6]) = -1.0;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[7]) = 1.0/(k*k*k*k*k*k*k)*(eta*k_*6.0+W1*(k*k*k)*3.0+W2*(k*k*k*k)*2.0+W2*(k*k*k*k*k*k*k))*4.0;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[8]) = 1.0;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[8]) = a;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[1]) = 1.0;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[3]) = -1.0;
  system.dC_dy.coeffRef(global_eqn_ids[4], global_var_ids[4]) = v*M_PI*(Ro+r)*8.0;
  system.dC_dy.coeffRef(global_eqn_ids[4], global_var_ids[5]) = M_PI*pow(Ro+r,2.0)*4.0;
  system.dC_dy.coeffRef(global_eqn_ids[5], global_var_ids[4]) = -1.0/(Ro*Ro)*(Ro*2.0+r*2.0);
  system.dC_dy.coeffRef(global_eqn_ids[5], global_var_ids[7]) = 1.0;
  system.dC_dy.coeffRef(global_eqn_ids[6], global_var_ids[5]) = -1.0;
  system.dC_dy.coeffRef(global_eqn_ids[7], global_var_ids[0]) = 1.0;
  system.dC_dy.coeffRef(global_eqn_ids[7], global_var_ids[2]) = -1.0;

  system.dC_dydot.coeffRef(global_eqn_ids[0], global_var_ids[5]) = d*rho;
  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[7]) = eta*(1.0/(k*k*k*k*k*k)*2.0-1.0)*-2.0;
  system.dC_dydot.coeffRef(global_eqn_ids[2], global_var_ids[8]) = 1.0;
  system.dC_dydot.coeffRef(global_eqn_ids[3], global_var_ids[9]) = -1.0;
  system.dC_dydot.coeffRef(global_eqn_ids[4], global_var_ids[9]) = -1.0;
  system.dC_dydot.coeffRef(global_eqn_ids[6], global_var_ids[4]) = 1.0;

  std::cout<<"wtf2"<<std::endl;
}

//double BloodVesselNew::get_active_stress(std::vector<double> &parameters) {
//  return 5.0;
//}%%%%%%%%%%%%%%%%%% come back to add active stress

// void BloodVesselNew::update_gradient(
//     Eigen::SparseMatrix<double> &jacobian,
//     Eigen::Matrix<double, Eigen::Dynamic, 1> &residual,
//     Eigen::Matrix<double, Eigen::Dynamic, 1> &alpha, std::vector<double> &y,
//     std::vector<double> &dy) {
//   auto y0 = y[global_var_ids[0]];
//   auto y1 = y[global_var_ids[1]];
//   auto y2 = y[global_var_ids[2]];
//   auto y3 = y[global_var_ids[3]];
//   auto y4 = y[global_var_ids[4]];
//   auto y5 = y[global_var_ids[5]];
//   auto y6 = y[global_var_ids[6]];
//   auto y7 = y[global_var_ids[7]];
//   auto y8 = y[global_var_ids[8]];
//   auto y9 = y[global_var_ids[9]];
  
//   auto dy0 = dy[global_var_ids[0]];
//   auto dy1 = dy[global_var_ids[1]];
//   auto dy3 = dy[global_var_ids[3]];
//   auto dy4 = dy[global_var_ids[4]];
//   auto dy5 = dy[global_var_ids[5]];
//   auto dy6 = dy[global_var_ids[6]];
//   auto dy7 = dy[global_var_ids[7]];
//   auto dy8 = dy[global_var_ids[8]];
//   auto dy9 = dy[global_var_ids[9]];

//   auto rho = alpha[global_param_ids[ParamId::rho]];
//   auto d = alpha[global_param_ids[ParamId::d]];
//   auto Ro = alpha[global_param_ids[ParamId::Ro]];
//   auto W1 = alpha[global_param_ids[ParamId::W1]];
//   auto W2 = alpha[global_param_ids[ParamId::W2]];
//   auto eta = alpha[global_param_ids[ParamId::eta]];
//   auto a = alpha[global_param_ids[ParamId::a]];
//   auto sigma_o = alpha[global_param_ids[ParamId::sigma_o]];

  
//   jacobian.coeffRef(global_eqn_ids[0], global_param_ids[0]) = d*dy5;
//   jacobian.coeffRef(global_eqn_ids[0], global_param_ids[1]) = rho*dy5+(y6*(y4/Ro+1.0))/Ro;
//   jacobian.coeffRef(global_eqn_ids[0], global_param_ids[2]) = 1.0/(Ro*Ro*Ro)*(y2*(y4*y4)*2.0-Ro*y6*d+y2*Ro*y4*2.0-y6*d*y4*2.0);
//   jacobian.coeffRef(global_eqn_ids[1], global_param_ids[3]) = 1.0/(y7*y7*y7)*-4.0+4.0;
//   jacobian.coeffRef(global_eqn_ids[1], global_param_ids[4]) = 1.0/(y7*y7)*(y7*y7*y7-1.0)*4.0;
//   jacobian.coeffRef(global_eqn_ids[1], global_param_ids[5]) = dy7*(1.0/(y7*y7*y7*y7*y7*y7)*2.0-1.0)*-2.0;
//   jacobian.coeffRef(global_eqn_ids[2], global_param_ids[6]) = -sigma_o+y8;
//   jacobian.coeffRef(global_eqn_ids[2], global_param_ids[7]) = -a;
//   jacobian.coeffRef(global_eqn_ids[4], global_param_ids[2]) = y5*M_PI*(Ro+y4)*8.0;
//   jacobian.coeffRef(global_eqn_ids[5], global_param_ids[2]) = 1.0/(Ro*Ro*Ro)*y4*(Ro+y4)*2.0;

//   residual(global_eqn_ids[0]) =
//       -y2 * pow(y4/Ro+1.0,2.0) + d * rho * dy5 + (y6 * d * (y4/Ro + 1.0))/Ro;
//   residual(global_eqn_ids[1]) =
//       -y6 + y8 - (1.0/(y7*y7*y7) * 4.0 - 4.0) * (W1 + W2 * y7) - eta * dy7 * (1.0/(y7*y7*y7*y7*y7*y7) * 2.0 - 1.0) * 2.0;
//   residual(global_eqn_ids[2]) =
//       dy8 - a * sigma_o + a * y8;
//   residual(global_eqn_ids[3]) =
//       y1 - y3 - dy9;
//   residual(global_eqn_ids[4]) =
//       -dy9 + (Ro * Ro) * y5 * M_PI * pow(y4/Ro + 1.0,2.0) * 4.0;
//   residual(global_eqn_ids[5]) =
//       y7 - pow(y4/Ro + 1.0,2.0);
//   residual(global_eqn_ids[6]) =
//       dy4 - y5;
//   residual(global_eqn_ids[7]) =
//       y0 - y2;
// }
