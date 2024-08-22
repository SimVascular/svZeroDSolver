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

#include "OpenLoopCoronaryBC.h"

void OpenLoopCoronaryBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {"volume_im"});
}

void OpenLoopCoronaryBC::update_constant(SparseSystem &system,
                                         std::vector<double> &parameters) {
  auto Ra = parameters[global_param_ids[0]];
  auto Ram = parameters[global_param_ids[1]];
  auto Rv = parameters[global_param_ids[2]];
  auto Ca = parameters[global_param_ids[3]];
  auto Cim = parameters[global_param_ids[4]];

  if (steady) {
    // Different equations for steady initial condition
    // Equations are:
    // -P_in + (Ra+Ram+Rv)Q_in + Pv = 0
    // V_im = 0
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = -1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = Ra + Ram + Rv;
  } else {
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = Cim * Rv;
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = Cim * Rv;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = -Cim * Rv * Ra;
    system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -(Rv + Ram);

    system.E.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -Ca * Cim * Rv;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
        Ra * Ca * Cim * Rv;
    system.E.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -Cim * Rv;
    system.E.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -Cim * Rv * Ram;
  }
}

void OpenLoopCoronaryBC::update_time(SparseSystem &system,
                                     std::vector<double> &parameters) {
  auto Ram = parameters[global_param_ids[1]];
  auto Rv = parameters[global_param_ids[2]];
  auto Cim = parameters[global_param_ids[4]];
  auto Pim = parameters[global_param_ids[5]];
  auto Pv = parameters[global_param_ids[6]];

  if (steady) {
    system.C(global_eqn_ids[1]) = Pv;
  } else {
    system.C(global_eqn_ids[0]) = Cim * (-Pim + Pv + this->Pim_0 - this->P_Cim_0);
    system.C(global_eqn_ids[1]) = (Ram * Cim * Pv) - Cim * (Rv + Ram) * (Pim + this->P_Cim_0 - this->Pim_0);
  }
}

void OpenLoopCoronaryBC::setup_initial_state_dependent_params(State initial_state, std::vector<double> &parameters) {
  auto P_in = initial_state.y[global_var_ids[0]];
  auto Q_in = initial_state.y[global_var_ids[1]];
  auto P_in_dot = initial_state.ydot[global_var_ids[0]];
  auto Q_in_dot = initial_state.ydot[global_var_ids[1]];
  auto Ra = parameters[global_param_ids[0]];
  auto Ram = parameters[global_param_ids[1]];
  auto Ca = parameters[global_param_ids[3]];
  // Pressure proximal to Ca and distal to Ra
  auto P_Ca = P_in - Ra*Q_in;
  auto P_Ca_dot = P_in_dot - Ra*Q_in_dot;
  // Flow into Ram (inflow minus flow into Ca)
  auto Q_am = Q_in - Ca*P_Ca_dot;
  // Pressure proximal to Cim/Vim and distal to Ram
  this->P_Cim_0 = P_Ca - Ram*Q_am;
  // Initial intramyocardial pressure
  this->Pim_0 = parameters[global_param_ids[5]];
}
