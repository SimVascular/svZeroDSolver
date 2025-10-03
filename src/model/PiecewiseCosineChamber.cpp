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

#include "PiecewiseCosineChamber.h"

void PiecewiseCosineChamber::setup_dofs(DOFHandler &dofhandler) {
  // Internal variable is chamber volume
  Block::setup_dofs_(dofhandler, 3, {"Vc"});
}

void PiecewiseCosineChamber::update_constant(SparseSystem &system,
                                       std::vector<double> &parameters) {
  // Eq 0: P_in - E(t)(Vc - Vrest) = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;

  // Eq 1: P_in - P_out - L*dQ_out = 0
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -1.0;

  // Eq 2: Q_in - Q_out - dVc = 0
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[3]) = -1.0;
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[4]) = -1.0;
}

void PiecewiseCosineChamber::update_time(SparseSystem &system,
                                   std::vector<double> &parameters) {
  get_elastance_values(parameters);

  // Eq 0: P_in - E(t)(Vc - Vrest) = P_in - E(t)*Vc + E(t)*Vrest = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -1 * Elas;
  system.C.coeffRef(global_eqn_ids[0]) =
      Elas * parameters[global_param_ids[ParamId::VREST]];
}

void PiecewiseCosineChamber::get_elastance_values(std::vector<double> &parameters) {
  double Emax = parameters[global_param_ids[ParamId::EMAX]];
  double Epass = parameters[global_param_ids[ParamId::EPASS]];
  double Vrest = parameters[global_param_ids[ParamId::VREST]];
  double contract_start = parameters[global_param_ids[ParamId::CSTART]];
  double relax_start = parameters[global_param_ids[ParamId::RSTART]];
  double contract_duration = parameters[global_param_ids[ParamId::CDUR]];
  double relax_duration = parameters[global_param_ids[ParamId::RDUR]];

  auto T_HB = model->cardiac_cycle_period;

  double phi = 0;

  auto piecewise_condition = fmod(model->time - contract_start, T_HB);

  if (0 <= piecewise_condition && piecewise_condition < contract_duration) {
    phi = 0.5 * (1 - cos((M_PI * piecewise_condition) / contract_duration));
  } else {
    piecewise_condition = fmod(model->time - relax_start, T_HB);
    if (0 <= piecewise_condition && piecewise_condition < relax_duration) {
      phi = 0.5 * (1 + cos((M_PI * piecewise_condition) / relax_duration));
    }
  }

  Elas = Epass + Emax * phi;
}
