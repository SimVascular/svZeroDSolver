// This file is part of svZeroDSolver licensed under the MIT License
// See the LICENSE.md file for license information

#include "ClosedLoopCoronaryLeftBC.h"

#include "Model.h"

void ClosedLoopCoronaryLeftBC::setup_model_dependent_params() {
  auto heart_block = model->get_block("CLH");
  im_param_id =
      heart_block->global_param_ids[ClosedLoopHeartPulmonary::ParamId::IML];
  ventricle_var_id =
      heart_block->global_var_ids[13];  // Solution ID for LV pressure
}
