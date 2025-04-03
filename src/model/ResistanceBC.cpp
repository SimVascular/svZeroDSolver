// This file is part of svZeroDSolver licensed under the MIT License
// See the LICENSE.md file for license information

#include "ResistanceBC.h"

void ResistanceBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 1, {});
}

void ResistanceBC::update_constant(SparseSystem &system,
                                   std::vector<double> &parameters) {
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
}

void ResistanceBC::update_time(SparseSystem &system,
                               std::vector<double> &parameters) {
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
      -parameters[global_param_ids[0]];
  system.C(global_eqn_ids[0]) = -parameters[global_param_ids[1]];
}
