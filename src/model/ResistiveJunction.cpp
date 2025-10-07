// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ResistiveJunction.h"

void ResistiveJunction::setup_dofs(DOFHandler &dofhandler) {
  // Set number of equations of a junction block based on number of
  // inlets/outlets. Must be set before calling parent constructor
  num_inlets = inlet_nodes.size();
  num_outlets = outlet_nodes.size();
  Block::setup_dofs_(dofhandler, num_inlets + num_outlets + 1, {"pressure_c"});
  num_triplets.F = (num_inlets + num_outlets) * 4;
}

void ResistiveJunction::update_constant(SparseSystem &system,
                                        std::vector<double> &parameters) {
  for (size_t i = 0; i < num_inlets; i++) {
    system.F.coeffRef(global_eqn_ids[i], global_var_ids[i * 2]) = 1.0;
    system.F.coeffRef(global_eqn_ids[i], global_var_ids[i * 2 + 1]) =
        -parameters[global_param_ids[i]];
    system.F.coeffRef(global_eqn_ids[i], global_var_ids.back()) = -1.0;
  }

  for (size_t i = num_inlets; i < num_inlets + num_outlets; i++) {
    system.F.coeffRef(global_eqn_ids[i], global_var_ids[i * 2]) = -1.0;
    system.F.coeffRef(global_eqn_ids[i], global_var_ids[i * 2 + 1]) =
        -parameters[global_param_ids[i]];
    system.F.coeffRef(global_eqn_ids[i], global_var_ids.back()) = 1.0;
  }

  // Mass conservation
  for (size_t i = 1; i < num_inlets * 2; i = i + 2) {
    system.F.coeffRef(global_eqn_ids[num_inlets + num_outlets],
                      global_var_ids[i]) = 1.0;
  }

  for (size_t i = (num_inlets * 2) + 1; i < (num_inlets + num_outlets) * 2;
       i = i + 2) {
    system.F.coeffRef(global_eqn_ids[num_inlets + num_outlets],
                      global_var_ids[i]) = -1.0;
  }
}
