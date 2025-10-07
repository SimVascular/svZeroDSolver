// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
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
