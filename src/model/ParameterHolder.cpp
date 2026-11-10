// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ParameterHolder.h"

ParameterHolder::ParameterHolder(
    const std::vector<std::pair<std::string, InputParameter>>& props)
    : input_param_properties(props) {
  for (const auto& p : props) {
    if (p.second.is_number) {
      params_[p.first] = p.second.is_optional ? p.second.default_val : 0.0;
    }
  }
}

void ParameterHolder::set_param(const std::string& name, double value) {
  params_[name] = value;
}

void ParameterHolder::set_param_vector(const std::string& name,
                                       const std::vector<double>& value) {
  params_vec_[name] = value;
}
