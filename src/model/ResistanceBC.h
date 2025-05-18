// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ResistanceBC.h
 * @brief model::ResistanceBC source file
 */
#ifndef SVZERODSOLVER_MODEL_RESISTANCEBC_HPP_
#define SVZERODSOLVER_MODEL_RESISTANCEBC_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief Resistance boundary condition.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1.0,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R$, *-*] (3,0)
 * node[anchor=south]{$P_{d}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_{in}-P_d=R \cdot Q_{in}
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{ll}P_{in} & Q_{in}\end{array}\right]^{T}
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ll}1 & -R\end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{C}^{e}=\left[-P_d\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Resistance
 * * `1` Distal pressure
 *
 * ### Internal variables
 *
 * This block has no internal variables.
 *
 */
class ResistanceBC : public Block {
 public:
  /**
   * @brief Construct a new ResistanceBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ResistanceBC(int id, Model *model)
      : Block(id, model, BlockType::resistance_bc,
              BlockClass::boundary_condition,
              {{"R", InputParameter()}, {"Pd", InputParameter()}}) {}

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
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(SparseSystem &system, std::vector<double> &parameters);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{1, 0, 0};
};

#endif  // SVZERODSOLVER_MODEL_RESISTANCEBC_HPP_
