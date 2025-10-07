// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "WindkesselBC.h"

void WindkesselBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {"pressure_c"});
}

void WindkesselBC::update_constant(SparseSystem &system,

                                   std::vector<double> &parameters) {
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = -1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -1.0;
}

void WindkesselBC::update_time(SparseSystem &system,

                               std::vector<double> &parameters) {
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[2]) =
      -parameters[global_param_ids[2]] * parameters[global_param_ids[1]];
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
      -parameters[global_param_ids[0]];
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) =
      parameters[global_param_ids[2]];
  system.C(global_eqn_ids[1]) = parameters[global_param_ids[3]];
}
