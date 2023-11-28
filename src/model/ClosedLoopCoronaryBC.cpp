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

#include "ClosedLoopCoronaryBC.h"

#include "Model.h"

void ClosedLoopCoronaryBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 3, {"volume_im"});
}

void ClosedLoopCoronaryBC::update_constant(SparseSystem &system,
                                           std::vector<double> &parameters) {
  auto ra = parameters[global_param_ids[ParamId::RA]];
  auto ram = parameters[global_param_ids[ParamId::RAM]];
  auto rv = parameters[global_param_ids[ParamId::RV]];
  auto ca = parameters[global_param_ids[ParamId::CA]];
  auto cim = parameters[global_param_ids[ParamId::CIM]];

  system.E.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -ram * ca;
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[1]) = ram * ra * ca;
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[0]) = -ca;
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[1]) = ca * ra;
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[4]) = -1.0;

  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = -1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = (ra + ram);
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[3]) = rv;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[2]) = cim;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[3]) = cim * rv;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[4]) = -1.0;
}

void ClosedLoopCoronaryBC::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {
  auto cim = parameters[global_param_ids[ParamId::CIM]];
  auto im = parameters[im_param_id];
  auto pim = im * y[ventricle_var_id];
  system.C(global_eqn_ids[2]) = -cim * pim;
}
