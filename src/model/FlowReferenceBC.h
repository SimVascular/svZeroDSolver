// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file FlowReferenceBC.h
 * @brief model::FlowReferenceBC source file
 */
#ifndef SVZERODSOLVER_MODEL_FLOWREFERENCEBC_HPP_
#define SVZERODSOLVER_MODEL_FLOWREFERENCEBC_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief Flow reference boundary condition.
 *
 * Applies a prescribed flow to a boundary.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$\hat{Q}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P$} to [short, *-] (1.2,0) ;
 * \draw [-latex] (1.4,0) -- (2.2,0) node[right] {$Q$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * Q=\hat{Q}
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{ll}P^{e} & Q^{e}\end{array}\right]^{T}
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ll}0 & 1\end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{C}^{e}=\left[\hat{Q}\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Flow
 *
 * ### Internal variables
 *
 * This block has no internal variables.
 *
 */
class FlowReferenceBC : public Block {
 public:
  /**
   * @brief Construct a new FlowReferenceBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  FlowReferenceBC(int id, Model *model)
      : Block(id, model, BlockType::flow_bc, BlockClass::boundary_condition,
              {{"t", InputParameter(false, true)},
               {"Q", InputParameter(false, true)}}) {}

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

#endif  // SVZERODSOLVER_MODEL_FLOWREFERENCEBC_HPP_
