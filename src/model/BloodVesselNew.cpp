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
  Block::setup_dofs_(dofhandler, 7, {"r", "v", "S", "tau", "V"});
}

void BloodVesselNew::update_constant(SparseSystem &system,
                                     std::vector<double> &parameters) {
  double rho = parameters[global_param_ids[ParamId::rho]];
  double d = parameters[global_param_ids[ParamId::d]];

  system.E.coeffRef(global_eqn_ids[0], global_var_ids[5]) = d * rho;
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[7]) = 1.0;

  //
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[6]) = -1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[7]) = 1.0;

  // Qin - Qout - dV/dt = 0
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[3]) = -1.0;
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[8]) = -1.0;

  system.F.coeffRef(global_eqn_ids[4], global_var_ids[8]) = -1.0;

  // dr/dt - v = 0
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[4]) = 1.0;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[5]) = -1.0;

  // Pin - Pout = 0
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[2]) = -1.0;
}

void BloodVesselNew::update_time(SparseSystem &system,
                                 std::vector<double> &parameters) {
  get_elastance_values(parameters);
  system.C(global_eqn_ids[2]) = -a_plus * sigma_o;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[7]) = a;
}

void BloodVesselNew::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {
  double d = parameters[global_param_ids[ParamId::d]];
  double Ro = parameters[global_param_ids[ParamId::Ro]];
  double W1 = parameters[global_param_ids[ParamId::W1]];
  double W2 = parameters[global_param_ids[ParamId::W2]];
  double eta = parameters[global_param_ids[ParamId::eta]];
  double sigma_o = parameters[global_param_ids[ParamId::sigma_o]];

  // y = [Pin, Qin, Pout, Qout, r, v, S, tau, V]
  //     [0,     1,    2,    3, 4, 5, 6,   7, 8]

  double Pout = y[global_var_ids[2]];
  double r = y[global_var_ids[4]];
  double v = y[global_var_ids[5]];
  double S = y[global_var_ids[6]];
  double tau = y[global_var_ids[7]];

  double r_ = dy[global_var_ids[4]];
  double V_ = dy[global_var_ids[8]];

  double rRo = r / Ro;
  double CCsqr = rRo + 1.0;
  double CC = CCsqr * CCsqr;
  double W_term = W1 + W2 * CC;
  double power_term = 4.0 * (1.0 / pow(CC, 3) - 1.0);
  double eta_term =
      (2.0 / pow(CC, 6) - 1.0) * (2.0 * CCsqr) * 2.0 * eta * r_ / Ro;

  system.C(global_eqn_ids[0]) = -Pout * CC + (S * d * CCsqr) / Ro;
  system.C(global_eqn_ids[1]) = -power_term * W_term - eta_term;
  system.C(global_eqn_ids[4]) = -V_ + 4.0 * M_PI * Ro * Ro * v * CC;

  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -CC;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      -2.0 * (Pout * CCsqr - S * d / 2.0) / Ro;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[6]) = d * CCsqr / Ro;

  double common_factor = 1.0 / pow(CCsqr, 7);
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4]) =
      24.0 * Ro * Ro * common_factor * (W1 + W2 * CC) +
      96.0 * pow(Ro, 10) * eta * r_ / pow(CC, 6) -
      8.0 * W2 * CCsqr * (pow(CC, -3) - 1.0) -
      4.0 * eta * r_ * (2.0 / pow(CC, 6) - 1.0);
  system.dC_dy.coeffRef(global_eqn_ids[4], global_var_ids[4]) =
      8.0 * M_PI * v * (Ro + r);
  system.dC_dy.coeffRef(global_eqn_ids[4], global_var_ids[5]) = 4.0 * M_PI * CC;
  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[4]) =
      -4.0 * eta * CCsqr * (2.0 / pow(CC, 6) - 1.0);
}

void BloodVesselNew::get_elastance_values(std::vector<double> &parameters) {
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

  double f = S_plus * S_minus;
  double a_t = alpha_max * f + alpha_min * (1 - f);

  a = std::abs(a_t);
  a_plus = std::max(a_t, 0.0);
}
