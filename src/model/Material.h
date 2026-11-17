// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file Material.h
 * @brief Hyperelastic wall materials shared by the chamber blocks
 */

#ifndef SVZERODSOLVER_MODEL_MATERIAL_HPP_
#define SVZERODSOLVER_MODEL_MATERIAL_HPP_

#include <Eigen/Core>
#include <complex>
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
 * The stress is available for real and complex arguments, so that blocks can
 * compute tangents by complex step. The invariants are
 * \f$\bar I_1 = J^{-2/3}\,\mathrm{tr}\,C\f$, \f$\bar I_2 = J^{-4/3}\,
 * \tfrac12(\mathrm{tr}^2 C - \mathrm{tr}\,C^2)\f$ and
 * \f$\bar I_4 = J^{-2/3}\,\mathbf{e}_F\cdot C\cdot\mathbf{e}_F\f$.
 *
 * Types (JSON `material` object): `mooney_rivlin`, `exponential`,
 * `mooney_rivlin_exponential`.
 */
class Material : public ParameterHolder {
 public:
  using cplx = std::complex<double>;         ///< complex-step scalar
  using Mat3 = Eigen::Matrix<double, 3, 3>;  ///< real 3x3 matrix
  using Mat3c = Eigen::Matrix<cplx, 3, 3>;   ///< complex 3x3 matrix
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
   * @brief Deviatoric 2nd Piola-Kirchhoff stress for complex arguments
   * (complex-step tangents)
   *
   * @param C Right Cauchy-Green tensor
   * @param J Volume ratio
   * @param f Unit fiber direction
   * @return Stress \f$2\,\partial W/\partial C\f$
   */
  virtual Mat3c stress(const Mat3c& C, cplx J, const Vec3& f) const = 0;

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
  Mat3c stress(const Mat3c& C, cplx J, const Vec3& f) const override;

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
