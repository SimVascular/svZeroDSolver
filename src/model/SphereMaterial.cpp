// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "SphereMaterial.h"

#include <cmath>
#include <stdexcept>

SphereMaterial::SphereMaterial(
    const std::vector<std::pair<std::string, InputParameter>>& props)
    : input_param_properties(props) {
  for (const auto& p : props) {
    if (p.second.is_number) {
      params_[p.first] = p.second.is_optional ? p.second.default_val : 0.0;
    }
  }
}

void SphereMaterial::set_param(const std::string& name, double value) {
  params_[name] = value;
}

std::unique_ptr<SphereMaterial> SphereMaterial::create(
    const std::string& type_str) {
  if (type_str == "mooney_rivlin") {
    return std::make_unique<MooneyRivlinMaterial>();
  }
  if (type_str == "exponential") {
    return std::make_unique<ExponentialMaterial>();
  }
  throw std::runtime_error(
      "Unknown material type '" + type_str +
      "'. Must be one of: mooney_rivlin, exponential");
}

SphericalStressResult MooneyRivlinMaterial::compute(double radius,
                                                     double radius0) const {
  const double W1 = params_.at("W1");
  const double W2 = params_.at("W2");

  SphericalStressResult res;
  res.C_val =
      2 *
      (2 * pow(radius + radius0, 5) *
           (-pow(radius0, 6) + pow(radius + radius0, 6)) *
           (W1 * pow(radius0, 2) + W2 * pow(radius + radius0, 2))) /
      (pow(radius0, 2) * pow(radius + radius0, 11));

  res.dC_dy_radius =
      24 * W1 * pow(radius0, 6) / pow(radius + radius0, 7) +
      8 * W2 * radius / pow(radius0, 2) +
      16 * W2 * pow(radius0, 4) / pow(radius + radius0, 5) +
      8 * W2 / radius0;

  return res;
}

SphericalStressResult ExponentialMaterial::compute(double radius,
                                                    double radius0) const {
  const double C0 = params_.at("C0");
  const double C1 = params_.at("C1");
  const double C2 = params_.at("C2");
  const double C3 = params_.at("C3");

  // S2, S6: differences of the deformed/reference radius raised to the
  // matching power. A: the I1,iso-like bracket driving the C0/C1 term.
  const double S2 = pow(radius + radius0, 2) - pow(radius0, 2);
  const double S6 = pow(radius + radius0, 6) - pow(radius0, 6);
  const double A = pow(radius0, 6) -
                    3 * pow(radius0, 2) * pow(radius + radius0, 4) +
                    2 * pow(radius + radius0, 6);

  // Each exponential is evaluated identically several times in the
  // unexpanded expression, so it is computed once here.
  const double E1 =
      exp(C1 * A * A / (pow(radius0, 4) * pow(radius + radius0, 8)));
  const double E2 = exp(C3 * S2 * S2 / pow(radius0, 4));

  SphericalStressResult res;
  res.C_val =
      2 *
      (4 * C0 * C1 * (radius + radius0) * S6 * A * E1 +
       2 * C2 * C3 * pow(radius + radius0, 11) * S2 * E2) /
      (pow(radius0, 2) * pow(radius + radius0, 11));

  res.dC_dy_radius =
      2 *
      (-32 * C0 * pow(C1, 2) * S6 * (A * A) *
           (A - 3 * pow(radius + radius0, 4) * S2) * E1 +
       8 * C2 * pow(C3, 2) * pow(radius + radius0, 20) * (S2 * S2) * E2 +
       2 * pow(radius0, 4) * pow(radius + radius0, 8) *
           (12 * C0 * C1 * pow(radius + radius0, 6) * A * E1 +
            24 * C0 * C1 * pow(radius + radius0, 4) * S2 * S6 * E1 +
            2 * C0 * C1 * S6 * A * E1 +
            2 * C2 * C3 * pow(radius + radius0, 12) * E2 +
            11 * C2 * C3 * pow(radius + radius0, 10) * S2 * E2) +
       11 * pow(radius0, 4) * pow(radius + radius0, 7) *
           (-4 * C0 * C1 * (radius + radius0) * S6 * A * E1 -
            2 * C2 * C3 * pow(radius + radius0, 11) * S2 * E2)) /
      (pow(radius0, 6) * pow(radius + radius0, 19));

  return res;
}
