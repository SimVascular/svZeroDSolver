// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ClosedLoopCoronaryLeftBC.h"

#include "Model.h"

void ClosedLoopCoronaryLeftBC::setup_model_dependent_params() {
  auto heart_block = model->get_block("CLH");
  im_param_id =
      heart_block->global_param_ids[ClosedLoopHeartPulmonary::ParamId::IML];
  ventricle_var_id =
      heart_block->global_var_ids[13];  // Solution ID for LV pressure
}
