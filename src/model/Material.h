// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file Material.h
 * @brief Hyperelastic wall materials shared by the chamber blocks
 */

#ifndef SVZERODSOLVER_MODEL_MATERIAL_HPP_
#define SVZERODSOLVER_MODEL_MATERIAL_HPP_

#include <Eigen/Core>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Parameter.h"
#include "ParameterHolder.h"

/**
 * @brief Hyperelastic wall material of the chamber blocks
 *
 * A material is an isochoric strain energy \f$W(\bar C)\f$ of the invariants
 * of \f$\bar C = J^{-2/3} C\f$ and the fiber direction \f$\mathbf{e}_F\f$. It
 * returns the deviatoric 2nd Piola-Kirchhoff stress
 * \f$\Sigma^d = 2\,\partial W/\partial C\f$ for a right Cauchy-Green tensor
 * \f$C\f$ given in an orthonormal basis, with the fiber direction in the same
 * basis. The volumetric response belongs to the block: \ref ChamberSphere is
 * incompressible by its kinematics, \ref ChamberCylinder enforces
 * incompressibility with an element pressure.
 *
 * The invariants are
 * \f$\bar I_1 = J^{-2/3}\,\mathrm{tr}\,C\f$, \f$\bar I_2 = J^{-4/3}\,
 * \tfrac12(\mathrm{tr}^2 C - \mathrm{tr}\,C^2)\f$ and
 * \f$\bar I_4 = J^{-2/3}\,\mathbf{e}_F\cdot C\cdot\mathbf{e}_F\f$.
 *
 * ### Linearization
 *
 * For the exact tangents of the blocks, stress_derivatives() also returns the
 * directional derivatives of the stress,
 * \f[
 * \mathrm{D}\Sigma^d[\delta C, \delta J] =
 * \frac{\partial\Sigma^d}{\partial C} : \delta C
 * + \frac{\partial\Sigma^d}{\partial J}\,\delta J,
 * \f]
 * along given directions \f$(\delta C_k, \delta J_k)\f$. \f$C\f$ and \f$J\f$
 * are separate arguments since the block's kinematics give \f$J\f$: a block
 * passes \f$\delta C = 2\,\partial E/\partial\xi\f$ and
 * \f$\delta J = \partial J/\partial\xi\f$ for each of its kinematic unknowns
 * \f$\xi\f$. MooneyRivlinExponentialMaterial lists the formulas; the
 * derivatives of \f$W\f$ in the invariants come from
 * `scripts/chamber_tangents.py material`.
 *
 * Types (JSON `material` object): `mooney_rivlin`, `exponential`,
 * `mooney_rivlin_exponential`.
 */
class Material : public ParameterHolder {
 public:
  using Mat3 = Eigen::Matrix<double, 3, 3>;  ///< real 3x3 matrix
  using Vec3 = Eigen::Matrix<double, 3, 1>;  ///< real 3-vector

  /**
   * @brief Construct a Material
   *
   * @param props Parameter name/spec pairs for this material type
   */
  explicit Material(
      const std::vector<std::pair<std::string, InputParameter>>& props);

  /**
   * @brief Deviatoric 2nd Piola-Kirchhoff stress
   *
   * @param C Right Cauchy-Green tensor
   * @param J Volume ratio \f$J = \sqrt{\det C}\f$ (given by the block's
   * kinematics)
   * @param f Unit fiber direction
   * @return Stress \f$2\,\partial W/\partial C\f$
   */
  virtual Mat3 stress(const Mat3& C, double J, const Vec3& f) const = 0;

  /**
   * @brief Stress and its directional derivatives
   *
   * Computes \f$\mathrm{D}\Sigma^d[\delta C_k, \delta J_k]\f$ for
   * \f$k = 0,\dots,n-1\f$ (see the class documentation).
   *
   * @param C Right Cauchy-Green tensor
   * @param J Volume ratio
   * @param f Unit fiber direction
   * @param n Number of directions
   * @param dC Directions of \f$C\f$ (n tensors)
   * @param dJ Directions of \f$J\f$ (n values)
   * @param dS Derivatives of the stress along the directions (output, n
   * tensors)
   * @return Stress \f$2\,\partial W/\partial C\f$
   */
  virtual Mat3 stress_derivatives(const Mat3& C, double J, const Vec3& f, int n,
                                  const Mat3* dC, const double* dJ,
                                  Mat3* dS) const = 0;

  /**
   * @brief Factory: create a material from a type string
   *
   * @param type_str One of: "mooney_rivlin", "exponential",
   * "mooney_rivlin_exponential"
   * @return Unique pointer to the created material
   */
  static std::unique_ptr<Material> create(const std::string& type_str);
};

/**
 * @brief Mooney-Rivlin material with exponential isotropic and fiber terms
 *
 * \f[
 * W = W_1(\bar I_1-3) + W_2(\bar I_2-3) + C_0 e^{C_1(\bar I_1-3)^2}
 *   + C_2 e^{C_3(\bar I_4-1)^2},
 * \f]
 * the sum of the `mooney_rivlin` and `exponential` energies. It is the passive
 * energy of \cite genet23, whose constants \f$C_1,\dots,C_6\f$ are
 * \f$W_1, W_2, C_0, C_1, C_2, C_3\f$ here.
 *
 * The stress is
 * \f[
 * \Sigma^d = 2\left(w_1\frac{\partial\bar I_1}{\partial C}
 * + w_2\frac{\partial\bar I_2}{\partial C}
 * + w_4\frac{\partial\bar I_4}{\partial C}\right),
 * \f]
 * with the gradient \f$w_i = \partial W/\partial\bar I_i\f$,
 * \f$w_1 = W_1 + 2C_0C_1x_1e^{C_1x_1^2}\f$, \f$w_2 = W_2\f$ and
 * \f$w_4 = 2C_2C_3x_4e^{C_3x_4^2}\f$, where \f$x_1 = \bar I_1 - 3\f$ and
 * \f$x_4 = \bar I_4 - 1\f$, and the derivatives of the invariants
 * - \f$\partial\bar I_1/\partial C =
 *   J^{-2/3}(\mathbf{I} - \tfrac13 I_1 C^{-1})\f$,
 * - \f$\partial\bar I_2/\partial C =
 *   J^{-4/3}(I_1\mathbf{I} - C - \tfrac23 I_2 C^{-1})\f$,
 * - \f$\partial\bar I_4/\partial C =
 *   J^{-2/3}(\mathbf{e}_F\otimes\mathbf{e}_F - \tfrac13 I_4 C^{-1})\f$,
 *
 * with \f$I_1 = \mathrm{tr}\,C\f$, \f$I_2 = \tfrac12(I_1^2 - C:C)\f$,
 * \f$I_4 = \mathbf{e}_F\cdot C\cdot\mathbf{e}_F\f$ and
 * \f$\partial J/\partial C = \tfrac12 J C^{-1}\f$.
 *
 * Linearization: the Hessian of \f$W\f$ in the invariants has the nonzero
 * entries \f$h_{11} = 2C_0C_1(1 + 2C_1x_1^2)e^{C_1x_1^2}\f$ and
 * \f$h_{44} = 2C_2C_3(1 + 2C_3x_4^2)e^{C_3x_4^2}\f$, so
 * \f[
 * \mathrm{D}\Sigma^d = 2\left(
 * h_{11}\,\mathrm{D}\bar I_1\,\frac{\partial\bar I_1}{\partial C}
 * + h_{44}\,\mathrm{D}\bar I_4\,\frac{\partial\bar I_4}{\partial C}
 * + \sum_{i=1,2,4} w_i\,\mathrm{D}\frac{\partial\bar I_i}{\partial C}\right).
 * \f]
 * The product rule on the expressions above gives \f$\mathrm{D}\bar I_i\f$
 * and \f$\mathrm{D}(\partial\bar I_i/\partial C)\f$ from
 * - \f$\mathrm{D}I_1 = \mathrm{tr}\,\delta C\f$,
 *   \f$\mathrm{D}I_2 = I_1\,\mathrm{tr}\,\delta C - C:\delta C\f$,
 *   \f$\mathrm{D}I_4 = \mathbf{e}_F\cdot\delta C\cdot\mathbf{e}_F\f$,
 * - \f$\mathrm{D}C^{-1} = -C^{-1}\,\delta C\,C^{-1}\f$,
 * - \f$\mathrm{D}J^{-2/3} = -\tfrac23 J^{-5/3}\,\delta J\f$,
 *   \f$\mathrm{D}J^{-4/3} = -\tfrac43 J^{-7/3}\,\delta J\f$.
 *
 * The exponents \f$C_1x_1^2\f$ and \f$C_3x_4^2\f$ are capped at 300 so that
 * non-physical trial states stay finite. Above the cap the exponential is
 * constant, and \f$h_{11}\f$ and \f$h_{44}\f$ lose their terms
 * \f$2C_1x_1^2\f$ and \f$2C_3x_4^2\f$.
 *
 * Parameters: `W1`, `W2`, `C0`, `C1`, `C2`, `C3` (optional, default 0)
 */
class MooneyRivlinExponentialMaterial : public Material {
 public:
  MooneyRivlinExponentialMaterial()
      : MooneyRivlinExponentialMaterial(
            {{"W1", InputParameter(true, false, true, 0.0)},
             {"W2", InputParameter(true, false, true, 0.0)},
             {"C0", InputParameter(true, false, true, 0.0)},
             {"C1", InputParameter(true, false, true, 0.0)},
             {"C2", InputParameter(true, false, true, 0.0)},
             {"C3", InputParameter(true, false, true, 0.0)}}) {}

  Mat3 stress(const Mat3& C, double J, const Vec3& f) const override;
  Mat3 stress_derivatives(const Mat3& C, double J, const Vec3& f, int n,
                          const Mat3* dC, const double* dJ,
                          Mat3* dS) const override;

 protected:
  /**
   * @brief Construct with a subset of the parameters (the others are 0)
   *
   * @param props Parameter name/spec pairs
   */
  explicit MooneyRivlinExponentialMaterial(
      const std::vector<std::pair<std::string, InputParameter>>& props)
      : Material(props) {}

 private:
  /// Parameter value, 0 for a parameter the material type does not have
  double get(const std::string& name) const;
};

/**
 * @brief Mooney-Rivlin material (reduces to neo-Hookean for W2 = 0)
 *
 * \f[
 * W = W_1(\bar I_1-3) + W_2(\bar I_2-3)
 * \f]
 *
 * Parameters: `W1`, `W2`
 */
class MooneyRivlinMaterial : public MooneyRivlinExponentialMaterial {
 public:
  MooneyRivlinMaterial()
      : MooneyRivlinExponentialMaterial(
            {{"W1", InputParameter()}, {"W2", InputParameter()}}) {}
};

/**
 * @brief Exponential material (\cite caruel13)
 *
 * \f[
 * W = C_0 e^{C_1(\bar I_1-3)^2} + C_2 e^{C_3(\bar I_4-1)^2}
 * \f]
 *
 * Parameters: `C0`, `C1`, `C2`, `C3`
 */
class ExponentialMaterial : public MooneyRivlinExponentialMaterial {
 public:
  ExponentialMaterial()
      : MooneyRivlinExponentialMaterial({{"C0", InputParameter()},
                                         {"C1", InputParameter()},
                                         {"C2", InputParameter()},
                                         {"C3", InputParameter()}}) {}
};

#endif  // SVZERODSOLVER_MODEL_MATERIAL_HPP_
