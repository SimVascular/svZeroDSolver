// This file is part of svZeroDSolver licensed under Stanford University, The Regents of the University of 
//                                                   California, and others.
// 
// See the LICENSE.md file for license information

#include "PressureReferenceBC.h"

void PressureReferenceBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 1, {});
}

void PressureReferenceBC::update_constant(SparseSystem &system,
                                          std::vector<double> &parameters) {
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
}

void PressureReferenceBC::update_time(SparseSystem &system,
                                      std::vector<double> &parameters) {
  system.C(global_eqn_ids[0]) = -parameters[global_param_ids[0]];
}
