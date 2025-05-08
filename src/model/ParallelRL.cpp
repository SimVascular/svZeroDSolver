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

#include "ParallelRL.h"

void ParallelRL::setup_dofs(DOFHandler &dofhandler) {
  // 4 equations and 2 internal variables (Q_R and Q_L)
  Block::setup_dofs_(dofhandler, 2, {});
}

void ParallelRL::update_constant(SparseSystem &system,
                                std::vector<double> &parameters) {
  double resistance = parameters[global_param_ids[ParamId::RESISTANCE]];
  double inductance = parameters[global_param_ids[ParamId::INDUCTANCE]];

  // unknowns: y = [P_in, Q_in, P_out, Q_out]

  // R*Q_in - (P_in - P_out) - L dQ_in/dt = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = resistance;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 1.0;
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -inductance;
  
  // Q_in - Q_out = 0
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
}