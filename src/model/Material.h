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
 * of \f$\bar C = J^{-2/3} C\f$, with the fiber and sheet directions
 * \f$\mathbf{e}_F\f$, \f$\mathbf{e}_S\f$. It returns the deviatoric 2nd
 * Piola-Kirchhoff stress \f$\Sigma^d = 2\,\partial W/\partial C\f$ for a right
 * Cauchy-Green tensor \f$C\f$ given in an orthonormal basis, with the
 * directions in the same basis. The volumetric response belongs to the block:
 * \ref ChamberSphere is incompressible by its kinematics, \ref ChamberCylinder
 * enforces its volumetric law with an element pressure.
 *
 * The stress is available for real and complex arguments, so that blocks can
 * compute tangents by complex step. The invariants are
 * \f$\bar I_1 = J^{-2/3}\,\mathrm{tr}\,C\f$, \f$\bar I_2 = J^{-4/3}\,
 * \tfrac12(\mathrm{tr}^2 C - \mathrm{tr}\,C^2)\f$,
 * \f$\bar I_4 = J^{-2/3}\,\mathbf{e}_F\cdot C\cdot\mathbf{e}_F\f$,
 * \f$\bar I_{4s} = J^{-2/3}\,\mathbf{e}_S\cdot C\cdot\mathbf{e}_S\f$ and
 * \f$\bar I_8 = J^{-2/3}\,\mathbf{e}_F\cdot C\cdot\mathbf{e}_S\f$.
 *
 * Types (JSON `material` object): `mooney_rivlin`, `exponential`, `genet23`,
 * `holzapfel_ogden`.
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
   * @param s Unit sheet direction
   * @return Stress \f$2\,\partial W/\partial C\f$
   */
  virtual Mat3 stress(const Mat3& C, double J, const Vec3& f,
                      const Vec3& s) const = 0;

  /**
   * @brief Deviatoric 2nd Piola-Kirchhoff stress for complex arguments
   * (complex-step tangents)
   *
   * @param C Right Cauchy-Green tensor
   * @param J Volume ratio
   * @param f Unit fiber direction
   * @param s Unit sheet direction
   * @return Stress \f$2\,\partial W/\partial C\f$
   */
  virtual Mat3c stress(const Mat3c& C, cplx J, const Vec3& f,
                       const Vec3& s) const = 0;

  /**
   * @brief Factory: create a material from a type string
   *
   * @param type_str One of: "mooney_rivlin", "exponential", "genet23",
   * "holzapfel_ogden"
   * @return Unique pointer to the created material
   */
  static std::unique_ptr<Material> create(const std::string& type_str);
};

/**
 * @brief Energy of \cite genet23 (Eq. 24)
 *
 * \f[
 * W = C_1(\bar I_1-3) + C_2(\bar I_2-3) + C_3 e^{C_4(\bar I_1-3)^2}
 *   + C_5 e^{C_6(\bar I_4-1)^2}
 * \f]
 *
 * With \f$C_3=0\f$ it contains the svMultiPhysics laws `neoHookean`
 * (\f$C_1=\mu/2\f$), `MR` (\f$C_1,C_2\f$) and the fiber term of `HGO` without
 * dispersion (\f$C_5=a_4/(2b_4)\f$, \f$C_6=b_4\f$).
 *
 * Parameters: `C1`, `C2`, `C3`, `C4`, `C5`, `C6` (optional, default 0)
 */
class Genet23Material : public Material {
 public:
  Genet23Material()
      : Material({{"C1", InputParameter(true, false, true, 0.0)},
                  {"C2", InputParameter(true, false, true, 0.0)},
                  {"C3", InputParameter(true, false, true, 0.0)},
                  {"C4", InputParameter(true, false, true, 0.0)},
                  {"C5", InputParameter(true, false, true, 0.0)},
                  {"C6", InputParameter(true, false, true, 0.0)}}) {}

  Mat3 stress(const Mat3& C, double J, const Vec3& f,
              const Vec3& s) const override;
  Mat3c stress(const Mat3c& C, cplx J, const Vec3& f,
               const Vec3& s) const override;
};

/**
 * @brief Mooney-Rivlin material (reduces to neo-Hookean for W2 = 0)
 *
 * \f[
 * W = W_1(\bar I_1-3) + W_2(\bar I_2-3),
 * \f]
 * the `genet23` energy with \f$C_1=W_1\f$, \f$C_2=W_2\f$.
 *
 * Parameters: `W1`, `W2`
 */
class MooneyRivlinMaterial : public Material {
 public:
  MooneyRivlinMaterial()
      : Material({{"W1", InputParameter()}, {"W2", InputParameter()}}) {}

  Mat3 stress(const Mat3& C, double J, const Vec3& f,
              const Vec3& s) const override;
  Mat3c stress(const Mat3c& C, cplx J, const Vec3& f,
               const Vec3& s) const override;
};

/**
 * @brief Exponential material (\cite caruel13)
 *
 * \f[
 * W = C_0 e^{C_1(\bar I_1-3)^2} + C_2 e^{C_3(\bar I_4-1)^2},
 * \f]
 * the `genet23` energy with \f$C_3,C_4,C_5,C_6 = C_0,C_1,C_2,C_3\f$.
 *
 * Parameters: `C0`, `C1`, `C2`, `C3`
 */
class ExponentialMaterial : public Material {
 public:
  ExponentialMaterial()
      : Material({{"C0", InputParameter()},
                  {"C1", InputParameter()},
                  {"C2", InputParameter()},
                  {"C3", InputParameter()}}) {}

  Mat3 stress(const Mat3& C, double J, const Vec3& f,
              const Vec3& s) const override;
  Mat3c stress(const Mat3c& C, cplx J, const Vec3& f,
               const Vec3& s) const override;
};

/**
 * @brief Holzapfel-Ogden material of svMultiPhysics (`HolzapfelOgden`)
 *
 * \f[
 * W=\frac{a}{2b}e^{b(\bar I_1-3)}
 * +\sum_{i=f,s}\frac{a_i}{2b_i}\chi_i\left(e^{b_i(\bar I_{4i}-1)^2}-1\right)
 * +\frac{a_{fs}}{2b_{fs}}\left(e^{b_{fs}\bar I_8^2}-1\right)
 * \f]
 * with \f$\bar I_{4f}=\bar I_4\f$ and the smoothed Heaviside
 * \f$\chi_i=1/(1+e^{-k(\bar I_{4i}-1)})\f$.
 *
 * Parameters (svMultiPhysics names): `a`, `b`, `a4f`, `b4f`, `a4s`, `b4s`,
 * `afs`, `bfs`, `k` (optional, default 0). A fiber or sheet term is off when
 * its `a4f`/`a4s` is 0, otherwise its `b4f`/`b4s` must be positive.
 */
class HolzapfelOgdenMaterial : public Material {
 public:
  HolzapfelOgdenMaterial()
      : Material({{"a", InputParameter(true, false, true, 0.0)},
                  {"b", InputParameter(true, false, true, 0.0)},
                  {"a4f", InputParameter(true, false, true, 0.0)},
                  {"b4f", InputParameter(true, false, true, 0.0)},
                  {"a4s", InputParameter(true, false, true, 0.0)},
                  {"b4s", InputParameter(true, false, true, 0.0)},
                  {"afs", InputParameter(true, false, true, 0.0)},
                  {"bfs", InputParameter(true, false, true, 0.0)},
                  {"k", InputParameter(true, false, true, 0.0)}}) {}

  void finalize() override;

  Mat3 stress(const Mat3& C, double J, const Vec3& f,
              const Vec3& s) const override;
  Mat3c stress(const Mat3c& C, cplx J, const Vec3& f,
               const Vec3& s) const override;
};

#endif  // SVZERODSOLVER_MODEL_MATERIAL_HPP_
