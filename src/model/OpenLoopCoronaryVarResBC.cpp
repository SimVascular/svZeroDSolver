// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "OpenLoopCoronaryVarResBC.h"

#include <cmath>

#include "Model.h"

double OpenLoopCoronaryVarResBC::get_Ram(std::vector<double>& parameters,
                                         double time) const {
  // Extract parameters for time-varying resistance
  auto Ram_min = parameters[global_param_ids[6]];
  auto Ram_max = parameters[global_param_ids[7]];
  auto T_vc = parameters[global_param_ids[8]];
  auto T_vr = parameters[global_param_ids[9]];

  // Compute e(t) based on phase in cardiac cycle
  double e_t;
  if (time <= T_vc) {
    // Contraction phase
    e_t = 0.5 * (1.0 - cos(M_PI * time / T_vc));
  } else if (time <= T_vc + T_vr) {
    // Relaxation phase
    e_t = 0.5 * (1.0 + cos(M_PI * (time - T_vc) / T_vr));
  } else {
    // Rest phase
    e_t = 0.0;
  }

  // Compute time-varying resistance
  double sqrt_Ram_min = sqrt(Ram_min);
  double sqrt_Ram_max = sqrt(Ram_max);
  double Ram_t = sqrt_Ram_min + (sqrt_Ram_max - sqrt_Ram_min) * e_t;
  return Ram_t * Ram_t;
}
