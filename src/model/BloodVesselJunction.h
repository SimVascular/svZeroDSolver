// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file BloodVesselJunction.h
 * @brief model::BloodVesselJunction source file
 */
#ifndef SVZERODSOLVER_MODEL_BLOODVESSELJUNCTION_HPP_
#define SVZERODSOLVER_MODEL_BLOODVESSELJUNCTION_HPP_

#include "Block.h"
#include "BloodVessel.h"
#include "SparseSystem.h"

/**
 * @brief Junction between blood vessels
 *
 * Models a junction with one inlet and arbitrary outlets using
 * modified blood vessel elements between each inlet and outlet pair.
 *
 * \f[
 * \begin{circuitikz}
 * \draw node[left] {$Q_\text{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0.1) node[anchor=south]{$P_\text{in}$};
 * \draw (1,0) to [short, *-] (2.5,0.75);
 * \draw (1,0) to [short, *-] (2.5,-0.75);
 * \draw (2.5,0.75) node[anchor=south]{} to [generic, l_=$BV_{1}$, -*]
 * (4.5,0.75); \draw (2.4,0.75) node[anchor=south]{}; \draw (4.6,0.75)
 * node[anchor=south] {$P_{out,1}$}; \draw (2.5,-0.75) node[anchor=south]{}
 to
 * [generic, l^=$BV_{2}$, -*] (4.5,-0.75); \draw (2.4,-0.75)
 * node[anchor=north]{}; \draw (4.6,-0.75) node[anchor=north]
 * {$P_{out,2}$}; \draw [-latex] (4.7,0.75) -- (5.5,0.75) node[right]
 * {$Q_{out,1}$}; \draw [-latex] (4.7,-0.75) -- (5.5,-0.75) node[right]
 * {$Q_{out,2}$}; \end{circuitikz}
 * \f]
 *
 * Each blood vessel is modelled as:
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_\text{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_\text{in}$}
 * to [R, l=$R$, *-] (3,0)
 * to [R, l=$S$, -] (5,0)
 * (5,0) to [L, l=$L$, -*] (7,0)
 * node[anchor=south]{$P_\text{out}$};
 * \draw [-latex] (7.2,0) -- (8,0) node[right] {$Q_\text{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * Q_\text{in}-\sum_{i}^{n_{outlets}} Q_{out, i} = 0
 * \f]
 *
 * \f[
 * P_\text{in}-P_{\text{out},i} - (R+S|Q_{\text{out},i}|) \cdot Q_{\text{out},i}
 * - L \dot{Q}_{\text{out},i} = 0 \quad \forall i \in n_{outlets} \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lllllll}P_\text{in} & Q_\text{in}
 * & P_{out, 1} & Q_{out, 1} &
 * \dots & P_{out, i} & Q_{out, i}\end{array}\right] \f]
 *
 * \f[
 * \mathbf{F}^{e} = \left[\begin{array}{ccccccc}
 * 0 & 1 & 0 & -1 & 0 & -1 & \dots \\
 * 1 & 0 & -1 & -R_1 & 0 & 0 & \\
 * 1 & 0 & 0 & 0 & -1 & -R_2 & \\
 * \vdots & & & & & & \ddots
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{E}^{e} = \left[\begin{array}{ccccccc}
 * 0 & 0 & 0 & 0 & 0 & 0 & \dots \\
 * 0 & 0 & 0 & -L_1 & 0 & 0 & \\
 * 0 & 0 & 0 & 0 & 0 & -L_2 & \\
 * \vdots & & & & & & \ddots
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e} = \left[\begin{array}{c}
 * 0 \\
 * - S_1 |Q_{\text{in},1}| Q_{\text{in},1} \\
 * - S_2 |Q_{\text{in},2}| Q_{\text{in},2} \\
 * \vdots
 * \end{array}\right] \f]
 *
 * ### Gradient
 *
 * Gradient of the equations with respect to the parameters:
 *
 * \f[
 * \mathbf{J}^{e} = \left[\begin{array}{lllllllll}
 * 0 & 0 & \dots & 0 & 0 & \dots & 0 & 0 &
 * \dots \\
 * - y_4 & 0 & \dots & - \dot y_4 & 0 & \dots & |y_4| y_4 & 0 & \dots \\
 * 0 & - y_6 & \dots & 0 & - \dot y_6 & \dots & 0 & |y_6| y_6 & \dots \\
 * 0 & 0 & \ddots & 0 & 0 & \ddots & 0 & 0 & \ddots \\
 * \end{array}\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `i` Poiseuille resistance for inner blood vessel `i`
 * * `i+num_outlets` Inductance for inner blood vessel `i`
 * * `i+2*num_outlets` Stenosis coefficient for inner blood vessel `i`
 *
 * ### Internal variables
 *
 * This block has no internal variables.
 *
 */
class BloodVesselJunction : public Block {
 public:
  /**
   * @brief Construct a new BloodVesselJunction object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  BloodVesselJunction(int id, Model *model)
      : Block(id, model, BlockType::blood_vessel_junction, BlockClass::junction,
              {{"R_poiseuille", InputParameter()},
               {"L", InputParameter()},
               {"stenosis_coefficient", InputParameter()}}) {
    input_params_list = true;
  }

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
  void setup_dofs(DOFHandler &dofhandler);

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem &system, std::vector<double> &parameters);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  virtual void update_solution(
      SparseSystem &system, std::vector<double> &parameters,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

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
  void update_gradient(Eigen::SparseMatrix<double> &jacobian,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &residual,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &alpha,
                       std::vector<double> &y, std::vector<double> &dy);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets;

 private:
  int num_outlets;
};

#endif  // SVZERODSOLVER_MODEL_BLOODVESSELJUNCTION_HPP_
