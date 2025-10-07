// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file Junction.h
 * @brief model::Junction source file
 */
#ifndef SVZERODSOLVER_MODEL_JUNCTION_HPP_
#define SVZERODSOLVER_MODEL_JUNCTION_HPP_

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Junction
 *
 * Models a junction with arbitrary inlets and outlets. Across all inlets and
 * outlets of the junction, mass is conserved and pressure is continuous.
 *
 * Inlets and outlets can be specified in two ways. Either using `inlet_vessels`
 * and `outlet_vessels` keys in the JSON file, with the corresponding lists
 * specifying vessel IDs, or using `inlet_blocks` and `outlet_blocks` keys, with
 * the corresponding lists specifying the names of blocks as strings.
 *
 *
 * \f[
 * \begin{circuitikz}
 * \draw node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$} to [short, *-*] (3.0,0);
 * \draw (3,0) node[anchor=south]{} to [short, -*] (4.5,1.0);
 * \draw (4.3,1.1) node[anchor=south] {$P_{out,1}$};
 * \draw (3,0) node[anchor=south]{} to [short, -*] (4.5,-1.0);
 * \draw (4.3,-1.1) node[anchor=north] {$P_{out,2}$};
 * \draw [-latex] (4.65,1.1) -- (5.25,1.5) node[right] {$Q_{out,1}$};
 * \draw [-latex] (4.65,-1.1) -- (5.25,-1.5) node[right] {$Q_{out,2}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * \sum_{i}^{n_{inlets}} Q_{in, i}=\sum_{j}^{n_{outlets}} Q_{out, j}
 * \f]
 *
 * \f[
 * P_{i}=P_{j} \quad \mathrm{with} \quad i \neq j
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{llllllllll}P_{in, 1}^{e} & Q_{in, 1}^{e} &
 * \dots & P_{in, i}^{e} & Q_{in, i}^{e} & P_{out, 1}^{e} & Q_{out, 1}^{e} &
 * \dots & P_{out, i}^{e} & Q_{out, i}^{e}\end{array}\right] \f]
 *
 * Mass conservation
 *
 * \f[
 * \mathbf{F}^{e}_1 = \left[\begin{array}{llllllllll}0 & 1 & 0 & 1 & \dots & 0 &
 * -1 & 0 & -1 & \dots\end{array}\right] \f]
 *
 * Due to the pressure continuity, we can write for all independent pressure
 * pairs: \f[ \mathbf{F}^{e}_{2,...,n} = \left[\begin{array}{lllll}\dots &
 * \underbrace{1}_{P_i} & \dots & \underbrace{1}_{P_j} & \dots\end{array}\right]
 * \quad \mathrm{with} \quad i \neq j \f]
 *
 * ### Internal variables
 *
 * This block has no internal variables.
 *
 */
class Junction : public Block {
 public:
  static const BlockType block_type;    ///< Type of this block
  static const BlockClass block_class;  ///< Class of this block
  static const std::vector<InputParameter>
      input_params;  ///< List of input parameter names

  /**
   * @brief Construct a new Junction object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  Junction(int id, Model *model)
      : Block(id, model, BlockType::junction, BlockClass::junction, {}) {}
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
  TripletsContributions num_triplets{0, 0, 0};

  /**
   * @brief Number of inlets to the block.
   *
   */
  int num_inlets;
  /**
   * @brief Number of outlets from the block.
   *
   */
  int num_outlets;
};

#endif
