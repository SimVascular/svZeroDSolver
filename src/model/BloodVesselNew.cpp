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
#include "Model.h"


void BloodVesselNew::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 7, {"r","v","S","tau","V"});
}


void BloodVesselNew::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {
    
  get_elastance_values(parameters);  
  // Get parameters
  
  double rho = parameters[global_param_ids[ParamId::rho]];
  double d = parameters[global_param_ids[ParamId::d]];
  double Ro = parameters[global_param_ids[ParamId::Ro]];
  double W1 = parameters[global_param_ids[ParamId::W1]];
  double W2 = parameters[global_param_ids[ParamId::W2]];
  double eta = parameters[global_param_ids[ParamId::eta]];
  double sigma_o = parameters[global_param_ids[ParamId::sigma_o]];

  double Pin = y[global_var_ids[0]]; 
  double Qin = y[global_var_ids[1]]; 
  double Pout = y[global_var_ids[2]]; 
  double Qout = y[global_var_ids[3]]; 
  double r = y[global_var_ids[4]]; 
  double v = y[global_var_ids[5]]; 
  double S = y[global_var_ids[6]]; 
  double tau = y[global_var_ids[7]]; 
  double V = y[global_var_ids[8]]; 


  double Pin_ = dy[global_var_ids[0]]; 
  double Qin_ = dy[global_var_ids[1]]; 
  double Pout_ = dy[global_var_ids[2]]; 
  double Qout_ = dy[global_var_ids[3]]; 
  double r_ = dy[global_var_ids[4]]; 
  double v_ = dy[global_var_ids[5]]; 
  double S_ = dy[global_var_ids[6]]; 
  double tau_ = dy[global_var_ids[7]]; 
  double V_ = dy[global_var_ids[8]];
  

  // Set element contributions
  system.C(global_eqn_ids[0]) = -Pout*pow(r/Ro+1.0,2.0)+d*rho*v_+(S*d*(r/Ro+1.0))/Ro;
  system.C(global_eqn_ids[1]) = -S+tau-(1.0/pow(r/Ro+1.0,6.0)*4.0-4.0)*(W1+W2*pow(r/Ro+1.0,2.0))-(eta*r_*(1.0/pow(r/Ro+1.0,1.2E+1)*2.0-1.0)*((r*2.0)/Ro+2.0)*2.0)/Ro;
  system.C(global_eqn_ids[2]) = tau_-a_plus*sigma_o+a*tau;
  system.C(global_eqn_ids[3]) = Qin-Qout-V_;
  system.C(global_eqn_ids[4]) = -V_+(Ro*Ro)*v*M_PI*pow(r/Ro+1.0,2.0)*4.0;
  system.C(global_eqn_ids[5]) = r_-v;
  system.C(global_eqn_ids[6]) = Pin-Pout;

  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -pow(r/Ro+1.0,2.0);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -1.0/(Ro*Ro)*(Pout*Ro*2.0-S*d+Pout*r*2.0);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[6]) = 1.0/(Ro*Ro)*d*(Ro+r);
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4]) = (Ro*Ro*Ro*Ro)*1.0/pow(Ro+r,7.0)*((Ro*Ro)*W1+(Ro*Ro)*W2+W2*(r*r)+Ro*W2*r*2.0)*2.4E+1+pow(Ro,1.0E+1)*eta*r_*1.0/pow(Ro+r,1.2E+1)*9.6E+1-1.0/(Ro*Ro)*W2*(Ro+r)*((Ro*Ro*Ro*Ro*Ro*Ro)*1.0/pow(Ro+r,6.0)*4.0-4.0)*2.0-1.0/(Ro*Ro)*eta*r_*(pow(Ro,1.2E+1)*1.0/pow(Ro+r,1.2E+1)*2.0-1.0)*4.0;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[6]) = -1.0;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[7]) = 1.0;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[7]) = a;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[1]) = 1.0;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[3]) = -1.0;
  system.dC_dy.coeffRef(global_eqn_ids[4], global_var_ids[4]) = v*M_PI*(Ro+r)*8.0;
  system.dC_dy.coeffRef(global_eqn_ids[4], global_var_ids[5]) = M_PI*pow(Ro+r,2.0)*4.0;
  system.dC_dy.coeffRef(global_eqn_ids[5], global_var_ids[5]) = -1.0;
  system.dC_dy.coeffRef(global_eqn_ids[6], global_var_ids[0]) = 1.0;
  system.dC_dy.coeffRef(global_eqn_ids[6], global_var_ids[2]) = -1.0;

  system.dC_dydot.coeffRef(global_eqn_ids[0], global_var_ids[5]) = d*rho;
  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[4]) = 1.0/(Ro*Ro)*eta*(Ro+r)*(pow(Ro,1.2E+1)*1.0/pow(Ro+r,1.2E+1)*2.0-1.0)*-4.0;
  system.dC_dydot.coeffRef(global_eqn_ids[2], global_var_ids[7]) = 1.0;
  system.dC_dydot.coeffRef(global_eqn_ids[3], global_var_ids[8]) = -1.0;
  system.dC_dydot.coeffRef(global_eqn_ids[4], global_var_ids[8]) = -1.0;
  system.dC_dydot.coeffRef(global_eqn_ids[5], global_var_ids[4]) = 1.0;

  // std::cout << "Pin: " << Pin << std::endl;
  // std::cout << "Qin: " << Qin << std::endl;
  // std::cout << "Pout: " << Pout << std::endl;
  // std::cout << "Qout: " << Qout << std::endl;
  
}


void BloodVesselNew::get_elastance_values(
    std::vector<double> &parameters) {
  double alpha_max = parameters[global_param_ids[ParamId::alpha_max]];
  double alpha_min = parameters[global_param_ids[ParamId::alpha_min]];
  double tsys = parameters[global_param_ids[ParamId::tsys]];
  double tdias = parameters[global_param_ids[ParamId::tdias]];
  double steepness = parameters[global_param_ids[ParamId::steepness]];

  double t = model->time;

  auto T_cardiac = model->cardiac_cycle_period;
  auto t_in_cycle = fmod(model->time, T_cardiac);

  double S_plus = 0.5 * (1.0 + tanh((t_in_cycle - tsys) / steepness));
  double S_minus = 0.5 * (1.0 - tanh((t_in_cycle - tdias) / steepness));

  // // Activation indicator function f(t)
  double f = S_plus*S_minus;

  double a_t = alpha_max * f + alpha_min * (1 - f);

  a = std::abs(a_t);
  a_plus = std::max(a_t, 0.0); // |a(t)|+ = max(a(t), 0)
  

  std::cout << f << std::endl;
  std::cout << "S_plus: " << S_plus << std::endl;
  std::cout << "S_minus: " << S_minus << std::endl;
  std::cout << "time: " << t_in_cycle << std::endl;
  std::cout << "a: " << a << std::endl;
  std::cout << "a_plus: " << a_plus << std::endl;
  
  // if (t_in_cycle <= tsys) {
  //   f = 0.0; // Before systole
  // } else if (t_in_cycle <= tdias) {
  //   f = 1.0; // Plateau during systole
  // } else {
  //   f = 0.0; // After diastole
  // }

  // double S_plus = 0.5 * (1.0 + tanh((t_in_cycle - tsys) / steepness));
  // double S_minus = 0.5 * (1.0 - tanh((t_in_cycle - tdias) / steepness));

  // // // Activation indicator function f(t)
  // double f = S_plus*S_minus;

}

