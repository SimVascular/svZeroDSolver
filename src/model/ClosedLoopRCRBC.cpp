// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ClosedLoopRCRBC.h"

void ClosedLoopRCRBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 3, {"P_c"});
}

void ClosedLoopRCRBC::update_constant(SparseSystem &system,
                                      std::vector<double> &parameters) {
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[3]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[4]) = -1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[2]) = -1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[4]) = 1.0;

  // Below values can be unsteady if needed (not currently implemented)
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      parameters[global_param_ids[ParamId::C]];
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) =
      -parameters[global_param_ids[ParamId::RP]];
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[3]) =
      -parameters[global_param_ids[ParamId::RD]];
}
