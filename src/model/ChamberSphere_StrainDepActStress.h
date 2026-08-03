// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ChamberSphere_StrainDepActStress.h
 * @brief model::ChamberSphere_StrainDepActStress source file
 */
#ifndef SVZERODSOLVER_MODEL_ChamberSphere_StrainDepActStress_HPP_
#define SVZERODSOLVER_MODEL_ChamberSphere_StrainDepActStress_HPP_

#include <math.h>

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Spherical heart chamber model with strain dependent active stress
 *
 * Models the mechanical behavior of a spherical heart chamber with active
 * contraction. For reference, see \cite caruel13 Equations (13a-b) for
 * continuum mechanics and Equations (13c-f) and (3) for the active
 * contraction model.
 *
 * ### Helper Functions
 *
 * Cauchy-Green deformation tensor and time derivative:
 * \f[
 * C = \left(1 + \frac{r}{r_0} \right)^2
 * \f]
 * \f[
 * \dot{C} = 2 \left(1 + \frac{r}{r_0} \right) \frac{\dot{r}}{r_0}
 * \f]
 *
 * ### Governing equations
 *
 * 1. Balance of linear momentum:
 * \f[
 * \rho d_0 \dot{v} + \frac{d_0}{r_0} \left(1 + \frac{r}{r_0} \right) S -
 P_\text{out} C = 0
 * \f]
 *
 * 2. Spherical stress:
 * \f[
 * -S + \tau + 4 (1 - C^{-3}) (W_1 + C W_2) + 2 \eta \dot{C}
 * (1 - 2 C^{-6}) = 0
 * \f]
 *
 * 3. Volume change:
 * \f[
 * 4 \pi r_0^2 Cv - \dot{V} = 0
 * \f]
 *
 * 4. Acceleration:
 * \f[
 * \dot{r} - v = 0
 * \f]
 *
 * 5. Conservation of mass:
 * \f[
 * Q_\text{in} - Q_\text{out} - \dot{V} = 0
 * \f]
 *
 * 6. Pressure equality:
 * \f[
 * P_\text{in} - P_\text{out} = 0
 * \f]
 *
 * 7. Active stress:
 * \f[
 * \dot{\tau}=E_s \frac{e_{1D}-e_c}{(1+2e_c)^2}
 * \f]
 * with $e_{1D}$ being the strain in fiber direction
 *
 * and with evolution equations
 * \f[
 * \dot{\omega} = \frac{1}{\alpha_r}(m_0-\omega)
 * \f]
 * \f[
 * \dot{e_c} =
 \frac{1}{\mu}\left(E_s\frac{(e_{1D}-e_c)(1+2e_{1D})}{(1+2e_c)^3}-\tau_c\right)
 * \f]
 * \f[
 * \dot{k}_c = -(|\bar(u)|_{+}+\omega|\bar(u)|_{-}+\alpha|\dot{e}_c|)k_c+n_0 k_0
 |\bar(u)|_{+}
 * \f]
 * \f[
 * \dot{\tau}_c = -(|\bar(u)|_{+}+|\bar(u)|_{-}+\alpha|\dot{e}_c|)\tau_c + n_0
 \sigma_0 |\bar(u)|_{+} + k_c \dot{e}_c
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block:
 *
 * * `rho` - Density \f$\rho\f$
 * * `thick0` - Wall thickness \f$d_0\f$
 * * `radius0` - Reference radius \f$r_0\f$
 * * `W1` - Material constant \f$W_1\f$
 * * `W2` - Material constant \f$W_2\f$
 * * `eta` - Viscosity parameter \f$\eta\f$
 * * `E_s` - Spring stiffness (in series to contractile unit) \f$E_s\f$
 * * `mu` - Damping coefficient (in parallel to contractile unit) \f$\mu\f$
 * * `alpha_r` - Relaxation time constant \f$\alpha_r\f$
 * * `alpha` - Activation time constant \f$\alpha\f$
 * * `k_0` - Linear spring stiffness in contractile unit \f$k_0\f$
 * * `sigma_0` - Maximum active stress \f$\sigma_0\f$
 * * 'm_0' - Relaxation strain-dependence \f$m_0\f$
 * * `n_0` - Activation strain-dependence \f$n_0\f$
 * * `u_plus` - Positive part of activation function \f$u_+\f$
 * * `u_minus` - Negative part of activation function \f$u_-\f$
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `radius` - Chamber radius \f$r\f$
 * * `velo` - Chamber velocity \f$\dot{r}\f$
 * * `stress` - Spherical stress \f$S\f$
 * * `tau` - Active stress \f$\tau\f$
 * * `volume` - Chamber volume \f$V\f$
 *
 */
class ChamberSphere_StrainDepActStress : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    rho = 0,
    thick0 = 1,
    radius0 = 2,
    C0 = 3,
    C1 = 4,
    C2 = 5,
    C3 = 6,
    eta = 7,
    E_s = 8,
    mu = 9,
    alpha_r = 10,
    alpha = 11,
    k_0 = 12,
    sigma_0 = 13,
  };

  /**
   * @brief Construct a new ChamberSphere object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ChamberSphere_StrainDepActStress(int id, Model* model)
      : Block(id, model, BlockType::chamber_sphere, BlockClass::vessel,
              {{"rho", InputParameter()},
               {"thick0", InputParameter()},
               {"radius0", InputParameter()},
               {"C0", InputParameter()},
               {"C1", InputParameter()},
               {"C2", InputParameter()},
               {"C3", InputParameter()},
               {"eta", InputParameter()},
               {"E_s", InputParameter()},
               {"mu", InputParameter()},
               {"alpha_r", InputParameter()},
               {"alpha", InputParameter()},
               {"k_0", InputParameter()},
               {"sigma_0", InputParameter()}}) {}

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set \ref global_var_ids and \ref global_eqn_ids of the element based on the
   * number of equations and the number of internal variables of the
   * element.
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler& dofhandler);

  /**
   * @brief Update the constant contributions of the element in a sparse
   system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem& system, std::vector<double>& parameters);

  // /**
  //  * @brief Update the time-dependent contributions of the element in a
  //  sparse
  //  * system
  //  *
  //  * @param system System to update contributions at
  //  * @param parameters Parameters of the model
  //  */
  // void update_time(SparseSystem &system, std::vector<double> &parameters);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(SparseSystem& system, std::vector<double>& parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy);

  /**
   * @brief Update the active stress functions which depend on time and strain
   *
   * @param parameters Parameters of the model
   * @param e_c Current strain
   */
  void get_active_stress_values(std::vector<double>& parameters,
                                const double e_c);

 private:
  double n_0 = 0.0;      // activation stain-dependence
  double m_0 = 0.0;      // relaxation strain-dependence
  double u_plus = 0.0;   // positive part of activation function
  double u_minus = 0.0;  // negative part of activation function

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{0, 0, 18};
};

#endif  // SVZERODSOLVER_MODEL_ChamberSphere_StrainDepActStress_HPP_
