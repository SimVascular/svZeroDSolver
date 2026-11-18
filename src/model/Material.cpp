// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "Material.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>

Material::Material(
    const std::vector<std::pair<std::string, InputParameter>>& props)
    : ParameterHolder(props) {}

std::unique_ptr<Material> Material::create(const std::string& type_str) {
  if (type_str == "mooney_rivlin") {
    return std::make_unique<MooneyRivlinMaterial>();
  }
  if (type_str == "exponential") {
    return std::make_unique<ExponentialMaterial>();
  }
  if (type_str == "mooney_rivlin_exponential") {
    return std::make_unique<MooneyRivlinExponentialMaterial>();
  }
  throw std::runtime_error("Unknown material type '" + type_str +
                           "'. Must be one of: mooney_rivlin, exponential, "
                           "mooney_rivlin_exponential");
}

double MooneyRivlinExponentialMaterial::get(const std::string& name) const {
  const auto it = params_.find(name);
  return it == params_.end() ? 0.0 : it->second;
}

Material::Mat3 MooneyRivlinExponentialMaterial::stress(const Mat3& C, double J,
                                                       const Vec3& f) const {
  return stress_derivatives(C, J, f, 0, nullptr, nullptr, nullptr);
}

Material::Mat3 MooneyRivlinExponentialMaterial::stress_derivatives(
    const Mat3& C, double J, const Vec3& f, int n, const Mat3* dC,
    const double* dJ, Mat3* dS) const {
  const double W1 = get("W1"), W2 = get("W2"), C0 = get("C0"), C1 = get("C1"),
               C2 = get("C2"), C3 = get("C3");
  const Mat3 I = Mat3::Identity();
  const Mat3 M = f * f.transpose();
  const Mat3 Cinv = C.inverse();
  const double I1 = C.trace();
  const double I2 = 0.5 * (I1 * I1 - (C.array() * C.array()).sum());
  const double I4 = f.dot(C * f);
  const double Jm23 = std::pow(J, -2.0 / 3.0);
  const double Jm43 = Jm23 * Jm23;

  // Derivatives of the isochoric invariants with respect to C
  const Mat3 dI1b = Jm23 * (I - (I1 / 3.0) * Cinv);
  const Mat3 dI2b = Jm43 * (I1 * I - C - (2.0 / 3.0) * I2 * Cinv);
  const Mat3 dI4b = Jm23 * (M - (I4 / 3.0) * Cinv);

  // Gradient and Hessian of W in the isochoric invariants
  // (scripts/chamber_tangents.py material). The exponents are capped so that
  // non-physical trial states stay finite; above the cap the exponential is
  // constant.
  const double x1 = Jm23 * I1 - 3.0, x4 = Jm23 * I4 - 1.0;
  const double a1 = C1 * x1 * x1, a4 = C3 * x4 * x4;
  const double e1 = std::exp(std::min(a1, 300.0));
  const double e4 = std::exp(std::min(a4, 300.0));
  const double w1 = W1 + 2.0 * C0 * C1 * x1 * e1;
  const double w2 = W2;
  const double w4 = 2.0 * C2 * C3 * x4 * e4;
  const double h11 = 2.0 * C0 * C1 * e1 * (a1 < 300.0 ? 1.0 + 2.0 * a1 : 1.0);
  const double h44 = 2.0 * C2 * C3 * e4 * (a4 < 300.0 ? 1.0 + 2.0 * a4 : 1.0);

  for (int k = 0; k < n; k++) {
    const Mat3& D = dC[k];
    const double dI1 = D.trace();
    const double dI2 = I1 * dI1 - (C.array() * D.array()).sum();
    const double dI4 = f.dot(D * f);
    const double dJm23 = -(2.0 / 3.0) * Jm23 / J * dJ[k];
    const double dJm43 = 2.0 * Jm23 * dJm23;
    const Mat3 dCinv = -Cinv * D * Cinv;
    // directional derivatives of the invariants and of dI/dC
    const double dx1 = dJm23 * I1 + Jm23 * dI1;
    const double dx4 = dJm23 * I4 + Jm23 * dI4;
    const Mat3 ddI1b = dJm23 * (I - (I1 / 3.0) * Cinv) +
                       Jm23 * (-(dI1 / 3.0) * Cinv - (I1 / 3.0) * dCinv);
    const Mat3 ddI2b =
        dJm43 * (I1 * I - C - (2.0 / 3.0) * I2 * Cinv) +
        Jm43 * (dI1 * I - D - (2.0 / 3.0) * (dI2 * Cinv + I2 * dCinv));
    const Mat3 ddI4b = dJm23 * (M - (I4 / 3.0) * Cinv) +
                       Jm23 * (-(dI4 / 3.0) * Cinv - (I4 / 3.0) * dCinv);
    dS[k] = 2.0 * (h11 * dx1 * dI1b + w1 * ddI1b + w2 * ddI2b +
                   h44 * dx4 * dI4b + w4 * ddI4b);
  }
  return 2.0 * (w1 * dI1b + w2 * dI2b + w4 * dI4b);
}
