// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ActiveStress.h"

#include <algorithm>
#include <cmath>

ActiveStress::ActiveStress(
    const std::vector<std::pair<std::string, InputParameter>>& props)
    : input_param_properties(props) {
  for (const auto& p : props) {
    if (p.second.is_number) {
      params_[p.first] = p.second.is_optional ? p.second.default_val : 0.0;
    }
  }
}

void ActiveStress::set_param(const std::string& name, double value) {
  params_[name] = value;
}

ActiveStressResult ElastanceActiveStress::compute(double f, double alpha_max,
                                                   double alpha_min) const {
  const double act_t = alpha_max * f + alpha_min * (1 - f);

  ActiveStressResult res;
  res.act = std::abs(act_t);
  res.act_plus = std::max(act_t, 0.0);
  return res;
}
