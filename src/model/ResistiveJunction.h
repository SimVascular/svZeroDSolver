// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ResistiveJunction.h
 * @brief model::ResistiveJunction source file
 */
#ifndef SVZERODSOLVER_MODEL_RESISTIVEJUNCTION_HPP_
#define SVZERODSOLVER_MODEL_RESISTIVEJUNCTION_HPP_

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief ResistiveJunction
 *
 * Models a junction with arbitrary resistive inlets and outlets. Across all
 * inlets and outlets of the junction, mass is conserved.
 *
 * \f[
 * \begin{circuitikz}
 * \draw [-latex] (0.25,1.4) node[left] {$Q_{in,1}$} -- (0.85,1.1);
 * \draw [-latex] (0.25,-1.4) node[left] {$Q_{in,1}$} -- (0.85,-1.1);
 * \draw (1,1.0) node[anchor=south]{$P_{in,1}$} to [R, , l=$R_{in,1}$, *-*]
 * (3.0,0) node[anchor=north] {$P_{C}$}; \draw (1,-1.0)
 * node[anchor=north]{$P_{in, 2}$} to [R, , l=$R_{in,2}$, *-*] (3.0,0); \draw
 * (3,0) node[anchor=south]{} to [R, l=$R_{out,1}$, -*] (5,1.0); \draw
 (4.3,1.1)
 * node[anchor=south] {$P_{out,1}$}; \draw (3,0) node[anchor=south]{} to [R,
 * l=$R_{out,2}$, -*] (5,-1.0); \draw (4.3,-1.1) node[anchor=north]
 * {$P_{out,2}$}; \draw [-latex] (5.15,1.1) -- (5.75,1.4) node[right]
 * {$Q_{out,1}$}; \draw [-latex] (5.15,-1.1) -- (5.75,-1.4) node[right]
 * {$Q_{out,2}$}; \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * \sum_{i}^{n_{inlets}} Q_{in, i}=\sum_{j}^{n_{outlets}} Q_{out, j}
 * \f]
 *
 * \f[
 * P_{in,i}-P_{C}=R_{in,i} \cdot Q_{in,i}\quad \forall i\in n_{inlets}
 * \f]
 * \f[
 * P_{C}-P_{out,j}=R_{out,j} \cdot Q_{out,j}\quad \forall j\in n_{outlets}
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lllllllllll}P_{in, 1}^{e} & Q_{in, 1}^{e}
 * & \dots & P_{in, i}^{e} & Q_{in, i}^{e} & P_{out, 1}^{e} & Q_{out, 1}^{e} &
 * \dots & P_{out, i}^{e} & Q_{out, i}^{e} & P_{C}\end{array}\right] \f]
 *
 * Mass conservation
 *
 * \f[
 * \mathbf{F}^{e}_1 = \left[\begin{array}{lllllllllll}0 & 1 & 0 & 1 & \dots & 0
 * & -1 & 0 & -1 & \dots & 0\end{array}\right] \f]
 *
 * \f[ \mathbf{F}^{e}_{2,...,n} = \left[\begin{array}{lllll}\dots &
 * \underbrace{1}_{P_{in,i}} & \underbrace{-R_{in,i}}_{Q_{in,i}} & \dots &
 * \underbrace{-1}_{P_{C}}\end{array}\right] \quad \mathrm{with} \quad \forall
 * i\in n_{inlets}  \f]
 *
 * \f[ \mathbf{F}^{e}_{2,...,n} = \left[\begin{array}{lllll}\dots &
 * \underbrace{-1}_{P_{out,j}} & \underbrace{-R_{out,j}}_{Q_{out,j}} & \dots &
 * \underbrace{1}_{P_{C}}\end{array}\right] \quad \mathrm{with} \quad \forall
 * j\in n_{oulets}  \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `i` Poiseuille resistance for inner blood vessel `i`
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `pressure_c`: Pressure at the center of the junction
 *
 */
class ResistiveJunction : public Block {
 public:
  /**
   * @brief Construct a new ResistiveJunction object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ResistiveJunction(int id, Model *model)
      : Block(id, model, BlockType::resistive_junction, BlockClass::junction,
              {{"R", InputParameter()}}) {}

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
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{0, 0, 0};

 private:
  int num_inlets;
  int num_outlets;
};

#endif  // SVZERODSOLVER_MODEL_RESISTIVEJUNCTION_HPP_
