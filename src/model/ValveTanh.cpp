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

#include "ValveTanh.h"

void ValveTanh::setup_dofs(DOFHandler &dofhandler) {
  // set_up_dofs args: dofhandler (passed in), num equations, list of internal
  // variable names (strings) 2 eqns, one for Pressure, one for Flow
  Block::setup_dofs_(dofhandler, 2, {});
}

// update_constant updates matrices E and F from E(y,t)*y_dot + F(y,t)*y +
// c(y,t) = 0 with terms that DO NOT DEPEND ON THE SOLUTION
void ValveTanh::update_constant(SparseSystem &system,
                                std::vector<double> &parameters) {
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

// update_solution updates matrices E and F from E(y,t)*y_dot + F(y,t)*y +
// c(y,t) = 0 with terms that DO DEPEND ON THE SOLUTION (will change with each
// time step)
void ValveTanh::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {
  // Get states
  double p_in = y[global_var_ids[0]];
  double p_out = y[global_var_ids[2]];
  double q_in = y[global_var_ids[1]];
  // Get parameters
  double Rmin = parameters[global_param_ids[ParamId::RMIN]];
  double Rmax = parameters[global_param_ids[ParamId::RMAX]];
  double steep = parameters[global_param_ids[ParamId::STEEPNESS]];

  // Nonlinear term
  system.C(global_eqn_ids[0]) =
      -0.5 * q_in * (Rmax - Rmin) * tanh(steep * (p_out - p_in));

  // Derivatives of non-linear term
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[0]) =
      0.5 * q_in * (Rmax - Rmin) * steep *
      (1.0 - tanh(steep * (p_out - p_in)) * tanh(steep * (p_out - p_in)));
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
      -0.5 * (Rmax - Rmin) * tanh(steep * (p_out - p_in));
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) =
      -0.5 * q_in * (Rmax - Rmin) * steep *
      (1.0 - tanh(steep * (p_out - p_in)) * tanh(steep * (p_out - p_in)));
}
