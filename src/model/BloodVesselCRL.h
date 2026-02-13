// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file BloodVesselCRL.h
 * @brief model::BloodVesselCRL source file
 */
#ifndef SVZERODSOLVER_MODEL_BLOODVESSELCRL_HPP_
#define SVZERODSOLVER_MODEL_BLOODVESSELCRL_HPP_

#include <math.h>

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Capacitor-resistor-inductor blood vessel with optional stenosis
 *
 * Models the mechanical behavior of a bloodvesselCRL with optional stenosis.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R$, *-] (3,0)
 * to [R, l=$S$, -] (5,0)
 * (5,0) to [L, l=$L$, -*] (7,0)
 * node[anchor=south]{$P_{out}$}
 * (1,0) to [C, l=$C$, -] (1,-1.5)
 * node[ground]{};
 * \draw [-latex] (7.2,0) -- (8,0) node[right] {$Q_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_\text{in}-P_\text{out} - (R + S|Q_\text{out}|) Q_\text{out}-L
 * \dot{Q}_\text{out}=0 \f]
 *
 * \f[
 * Q_\text{in}-Q_\text{out} - C \dot{P}_\text{in}=0 \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{llll}P_{i n} & Q_{in} &
 * P_{out} & Q_{out}\end{array}\right]^\text{T} \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{cccc}
 * 1 & 0 & -1 &  -R \\
 * 0 &  1 &  0 & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{cccc}
 *  0 &  0 & 0 & -L \\
 * -C & 0 & 0 &  0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e} = S|Q_\text{out}|
 * \left[\begin{array}{c}
 * -Q_\text{out} \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \left(\frac{\partial\mathbf{c}}{\partial\mathbf{y}}\right)^{e} =
 *  S \text{sgn} (Q_\text{in})
 * \left[\begin{array}{cccc}
 * 0 & -2Q_\text{out}        & 0 & 0 \\
 * 0 & 0 & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \left(\frac{\partial\mathbf{c}}{\partial\dot{\mathbf{y}}}\right)^{e} =
 *  S|Q_\text{out}|
 * \left[\begin{array}{cccc}
 * 0 &  0 & 0 & 0 \\
 * 0 & 0 & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * with the stenosis resistance \f$ S=K_{t} \frac{\rho}{2
 * A_{o}^{2}}\left(\frac{A_{o}}{A_{s}}-1\right)^{2} \f$.
 * \f$R\f$, \f$C\f$, and \f$L\f$ refer to
 * Poisieuille resistance, capacitance and inductance, respectively.
 *
 * ### Gradient
 *
 * Gradient of the equations with respect to the parameters:
 *
 * \f[
 * \mathbf{J}^{e} = \left[\begin{array}{cccc}
 * -y_3 & 0 & -\dot{y}_3 & -|y_3|y_3 \\
 * 0 & 0 & -\dot{y}_0 & 0 \\
 * \end{array}\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Poiseuille resistance
 * * `1` Capacitance
 * * `2` Inductance
 * * `3` Stenosis coefficient
 *
 */
class BloodVesselCRL : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    RESISTANCE = 0,
    CAPACITANCE = 1,
    INDUCTANCE = 2,
    STENOSIS_COEFFICIENT = 3,
  };

  /**
   * @brief Construct a new BloodVesselCRL object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  BloodVesselCRL(int id, Model* model)
      : Block(id, model, BlockType::blood_vessel_CRL, BlockClass::vessel,
              {{"R_poiseuille", InputParameter()},
               {"C", InputParameter(true)},
               {"L", InputParameter(true)},
               {"stenosis_coefficient", InputParameter(true)}}) {}

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
   * @brief Set the gradient of the block contributions with respect to the
   * parameters
   *
   * @param jacobian Jacobian with respect to the parameters
   * @param alpha Current parameter vector
   * @param residual Residual with respect to the parameters
   * @param y Current solution
   * @param dy Time-derivative of the current solution
   */
  void update_gradient(Eigen::SparseMatrix<double>& jacobian,
                       Eigen::Matrix<double, Eigen::Dynamic, 1>& residual,
                       Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha,
                       std::vector<double>& y, std::vector<double>& dy);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{5, 3, 2};
};

#endif  // SVZERODSOLVER_MODEL_BLOODVESSELCRL_HPP_
