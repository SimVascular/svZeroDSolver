// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ClosedLoopCoronaryBC.h
 * @brief model::ClosedLoopCoronaryBC source file
 */
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYBC_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYBC_HPP_

#include "Block.h"
#include "ClosedLoopHeartPulmonary.h"
#include "SparseSystem.h"

/**
 * @brief Closed loop coronary boundary condition which is connected to
 * other blocks on both sides and the intramyocardial pressure is
 * specified by the pressure in a heart block (not as a parameter).
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_a$, *-] (3,0)
 * to [R, l=$R_{am}$, -] (5,0)
 * to [R, l=$R_v$, *-*] (7,0)
 * node[anchor=south]{$P_{out}$}
 * (5,0) to [C, l=$C_{im} \;V_{im}$, -*] (5,-1.5)
 * node[left]{$P_{im}$}
 * (3,0) to [C, l=$C_a$, -*] (3,-1.5)
 * node[left]{$P_a$};
 * \draw [-latex] (7.2,0) -- (8.0,0) node[right] {$Q_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_{out} - P_{in} + (R_{am}+R_a)Q_{in} + R_v Q_{out} + R_{am} C_a
 * \frac{dP_a}{dt} - R_{am} C_a \frac{dP_{in}}{dt} + R_{am} R_a C_a
 * \frac{dQ_{in}}{dt} = 0 \f]
 *
 * \f[
 * Q_{in} - Q_{out} + C_a \frac{dP_a}{dt} - C_a \frac{dP_{in}}{dt} + C_a R_a
 * \frac{dQ_{in}}{dt} - \frac{dV_{im}}{dt} = 0 \f]
 *
 * \f[
 * C_{im} P_{out} + C_{im} R_v Q_{out} - C_{im} P_{im} - V_{im} = 0
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lllll}P_{in} & Q_{in} & P_{out} &
 * Q_{out} & V_{im}\end{array}\right]^{T}, \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{ccccc}
 * -R_{am} C_{a} & R_{am} R_{a} C_{a} & 0 & 0 & 0 \\
 * -C_{a} & R_{a} C_{a} & 0 & 0 & -1 \\
 * 0 & 0 & 0 & 0 & 0 \\
 * \end{array}\right] \f]
 *
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccccc}
 * -1 & R_{am} + R_{a} & 1 & R_v & 0 \\
 * 0 & 1 & 0 & -1 & 0 \\
 * 0 & 0 & C_{im} & C_{im} R_v & -1 \\
 * \end{array}\right] \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * C_{a} R_{am} \frac{d P_{a}}{d t} \\
 * C_{a}\frac{d P_{a}}{d t} \\
 * -C_{im} P_{im}
 * \end{array}\right] \f]
 *
 * Assume \f$P_a=0\f$.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Ra: Small artery resistance
 * * `1` Ram: Microvascular resistance
 * * `2` Rv: Venous resistance
 * * `3` Ca: Small artery capacitance
 * * `4` Cim: Intramyocardial capacitance
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `volume_im`: Intramyocardial volume
 *
 */
class ClosedLoopCoronaryBC : public Block {
 public:
  /**
   * @brief Construct a ClosedLoopCoronaryBC object.
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   * @param block_type The specific type of block (left or right)
   */
  ClosedLoopCoronaryBC(int id, Model *model, BlockType block_type)
      : Block(id, model, block_type, BlockClass::closed_loop,
              {{"Ra", InputParameter()},
               {"Ram", InputParameter()},
               {"Rv", InputParameter()},
               {"Ca", InputParameter()},
               {"Cim", InputParameter()}}) {}

  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId { RA = 0, RAM = 1, RV = 2, CA = 3, CIM = 4 };

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
   * @brief Setup parameters that depend on the model
   *
   */
  virtual void setup_model_dependent_params() = 0;

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
  void update_solution(SparseSystem &system, std::vector<double> &parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{9, 5, 0};

 protected:
  int ventricle_var_id;  ///< Variable index of either left or right ventricle
  int im_param_id;       ///< Index of parameter Im
};

#endif  // SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYBC_HPP_
