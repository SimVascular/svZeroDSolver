// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "Material.h"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace {

using cplx = std::complex<double>;
template <typename T>
using Mat3 = Eigen::Matrix<T, 3, 3>;
using Vec3 = Eigen::Matrix<double, 3, 1>;

/// Real part of a scalar, for comparisons on complex-step arguments.
inline double re(double x) { return x; }
inline double re(const cplx& x) { return x.real(); }

/// Upper clamp compatible with complex-step: x below m, the constant m above.
template <typename T>
inline T clamp_max(T x, double m) {
  return re(x) < m ? x : T(m);
}

/// Stress of the energy
/// W = W1 (I1b - 3) + W2 (I2b - 3) + C0 exp(C1 (I1b - 3)^2)
///   + C2 exp(C3 (I4b - 1)^2)
/// with the coefficients c = {W1, W2, C0, C1, C2, C3}.
template <typename T>
Mat3<T> stress_t(const Mat3<T>& C, const T& J, const Vec3& f,
                 const double c[6]) {
  const Eigen::Matrix<T, 3, 1> ef = f.cast<T>();
  const Mat3<T> M = ef * ef.transpose();
  const Mat3<T> I = Mat3<T>::Identity();
  const Mat3<T> Cinv = C.inverse();
  const T I1 = C.trace();
  const T I2 = 0.5 * (I1 * I1 - (C.array() * C.array()).sum());
  const T I4 = (C * M).trace();
  const T Jm23 = std::pow(J, -2.0 / 3.0);
  const T Jm43 = Jm23 * Jm23;
  const T I1b = Jm23 * I1;
  const T I4b = Jm23 * I4;

  // Derivatives of the isochoric invariants with respect to C
  const Mat3<T> dI1b = Jm23 * (I - (I1 / 3.0) * Cinv);
  const Mat3<T> dI2b = Jm43 * (I1 * I - C - (2.0 / 3.0) * I2 * Cinv);
  const Mat3<T> dI4b = Jm23 * (M - (I4 / 3.0) * Cinv);

  // dW/d(isochoric invariant); the exponents are capped so that non-physical
  // trial states stay finite.
  const T a1 = c[3] * (I1b - 3.0) * (I1b - 3.0);
  const T a4 = c[5] * (I4b - 1.0) * (I4b - 1.0);
  const T w1 =
      c[0] + 2.0 * c[2] * c[3] * (I1b - 3.0) * std::exp(clamp_max(a1, 300.0));
  const T w2 = T(c[1]);
  const T w4 = 2.0 * c[4] * c[5] * (I4b - 1.0) * std::exp(clamp_max(a4, 300.0));
  return 2.0 * (w1 * dI1b + w2 * dI2b + w4 * dI4b);
}

}  // namespace

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
  throw std::runtime_error(
      "Unknown material type '" + type_str +
      "'. Must be one of: mooney_rivlin, exponential, "
      "mooney_rivlin_exponential");
}

double MooneyRivlinExponentialMaterial::get(const std::string& name) const {
  const auto it = params_.find(name);
  return it == params_.end() ? 0.0 : it->second;
}

Material::Mat3 MooneyRivlinExponentialMaterial::stress(const Mat3& C, double J,
                                                       const Vec3& f) const {
  const double c[6] = {get("W1"), get("W2"), get("C0"),
                       get("C1"), get("C2"), get("C3")};
  return stress_t<double>(C, J, f, c);
}

Material::Mat3c MooneyRivlinExponentialMaterial::stress(const Mat3c& C,
                                                        cplx J,
                                                        const Vec3& f) const {
  const double c[6] = {get("W1"), get("W2"), get("C0"),
                       get("C1"), get("C2"), get("C3")};
  return stress_t<cplx>(C, J, f, c);
}
