// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file SphereMaterial.h
 * @brief Material models for the ChamberSphere block
 */

#ifndef SVZERODSOLVER_MODEL_SPHEREMATERIAL_HPP_
#define SVZERODSOLVER_MODEL_SPHEREMATERIAL_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Parameter.h"

/**
 * @brief Return type for SphereMaterial::compute
 *
 * Bundles the nonlinear residual contribution and its derivatives that
 * appear in the spherical-stress equation of ChamberSphere.
 */
struct SphericalStressResult {
  double C_val;            ///< Residual contribution to C vector
  double dC_dy_radius;     ///< Derivative w.r.t. radius
  double dC_dydot_radius;  ///< Derivative w.r.t. d(radius)/dt
};

/**
 * @brief Base class for spherical chamber wall materials
 *
 * Subclasses implement the material-specific elastic stress term (including
 * viscous damping) in the spherical stress equation:
 * \f[
 * -S + \tau + f_\text{material}(r, r_0, \dot{r}) = 0
 * \f]
 *
 * Each subclass owns all its material parameters, including the viscosity
 * \f$\eta\f$.
 */
class SphereMaterial {
 public:
  /**
   * @brief Properties of the input parameters for this material
   * [(name, InputParameter), ...]
   */
  const std::vector<std::pair<std::string, InputParameter>> input_param_properties;

  /**
   * @brief Construct a SphereMaterial
   *
   * @param props Parameter name/spec pairs for this material type
   */
  SphereMaterial(
      const std::vector<std::pair<std::string, InputParameter>>& props);

  virtual ~SphereMaterial() = default;

  /**
   * @brief Compute material stress contributions at the current state
   *
   * @param radius     Radius perturbation \f$r\f$
   * @param radius0    Reference radius \f$r_0\f$
   * @param dradius_dt Time derivative \f$\dot{r}\f$
   * @return SphericalStressResult
   */
  virtual SphericalStressResult compute(double radius, double radius0,
                                        double dradius_dt) const = 0;

  /**
   * @brief Factory: create a material from a type string
   *
   * @param type_str One of: "mooney_rivlin", "exponential"
   * @return Unique pointer to the created material
   */
  static std::unique_ptr<SphereMaterial> create(const std::string& type_str);

  /**
   * @brief Set a scalar parameter value by name
   *
   * @param name  Parameter name
   * @param value Parameter value
   */
  void set_param(const std::string& name, double value);

 protected:
  /**
   * @brief Map of parameter names to their values
   */
  std::map<std::string, double> params_;
};

/**
 * @brief Mooney–Rivlin (neo-Hookean) material
 *
 * Implements:
 * \f[
 * f = 4(1 - C^{-3})(W_1 + C W_2) + \eta \dot{C}(1 + 2 C^{-6})
 * \f]
 *
 * Parameters: `W1`, `W2`, `eta`
 */
class MooneyRivlinMaterial : public SphereMaterial {
 public:
  MooneyRivlinMaterial()
      : SphereMaterial({{"W1", InputParameter()},
                        {"W2", InputParameter()},
                        {"eta", InputParameter()}}) {}

  SphericalStressResult compute(double radius, double radius0,
                                double dradius_dt) const override;
};

/**
 * @brief Exponential material
 *
 * Implements the exponential strain-energy density:
 * \f[
 * W_1(r) = C_0 \exp(C_1 [I_{1,\text{iso}}]^2),\quad
 * W_4(r) = C_2 \exp(C_3 [I_{4,\text{iso}}]^2)
 * \f]
 *
 * Parameters: `C0`, `C1`, `C2`, `C3`, `eta`
 */
class ExponentialMaterial : public SphereMaterial {
 public:
  ExponentialMaterial()
      : SphereMaterial({{"C0", InputParameter()},
                        {"C1", InputParameter()},
                        {"C2", InputParameter()},
                        {"C3", InputParameter()},
                        {"eta", InputParameter()}}) {}

  SphericalStressResult compute(double radius, double radius0,
                                double dradius_dt) const override;
};

#endif  // SVZERODSOLVER_MODEL_SPHEREMATERIAL_HPP_
