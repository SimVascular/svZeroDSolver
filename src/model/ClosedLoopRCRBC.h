// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ClosedLoopRCRbc.h
 * @brief model::ClosedLoopRCRBC source file
 */
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPRCRBC_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPRCRBC_HPP_

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Closed-loop RCR boundary condition.
 *
 * Models the mechanical behavior of a Windkessel boundary condition that is
 * connected to other blocks on both sides.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_p$, *-] (3,0)
 * node[anchor=south]{$P_{c}$}
 * to [R, l=$R_d$, *-*] (5,0)
 * node[anchor=south]{$P_{out}$}
 * (3,0) to [C, l=$C$, *-] (3,-1.5)
 * node[ground]{};
 * \draw [-latex] (5.2,0) -- (6.0,0) node[right] {$Q_{out}$} ;
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * C \frac{d P_c}{dt} + Q_{out} -Q_{in} = 0
 * \f]
 *
 * \f[
 * P_{in}-P_{c}-R_{p} Q_{in}=0
 * \f]
 *
 * \f[
 * P_{c} - P_{out} - R_{d} Q_{out}=0
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^e=\left[\begin{array}{lllll}P_{in} & Q_{in} & P_{out} & Q_{out} &
 * P_{c}\end{array}\right]^{T} \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{ccccc}
 * 0 & 0 & 0 & 0 & C \\
 * 0 & 0 & 0 & 0 & 0 \\
 * 0 & 0 & 0 & 0 & 0 \\
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccccc}
 * 0 & -1 & 1 & 0 & 0 \\
 * 1 & -R_p & 0 & 0 & -1 \\
 * 0 & 0 & -1 & -R_d & +1 \\
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * 0 \\
 * 0 \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Proximal resistance
 * * `1` Capacitance
 * * `2` Distal resistance
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `P_c`: Pressure at the capacitor
 *
 */
class ClosedLoopRCRBC : public Block {
 public:
  /**
   * @brief Construct a new ClosedLoopRCRBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ClosedLoopRCRBC(int id, Model *model)
      : Block(id, model, BlockType::closed_loop_rcr_bc,
              BlockClass::boundary_condition,
              {{"Rp", InputParameter()},
               {"C", InputParameter()},
               {"Rd", InputParameter()},
               {"closed_loop_outlet", InputParameter(true, false, false)}}) {}

  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    RP = 0,
    C = 1,
    RD = 2,
  };

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
   * @brief Update the constant contributions of the element in a sparse
   system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem &system, std::vector<double> &parameters);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{8, 1, 0};
};

#endif  // SVZERODSOLVER_MODEL_CLOSEDLOOPRCRBCBC_HPP_
