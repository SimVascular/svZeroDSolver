// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file PiecewiseValve.h
 * @brief model::PiecewiseValve source file
 */
#ifndef SVZERODSOLVER_MODEL_PIECEWISE_VALVE_HPP_
#define SVZERODSOLVER_MODEL_PIECEWISE_VALVE_HPP_

#include <math.h>

#include "Block.h"
#include "SparseSystem.h"
#include "debug.h"

/**
 * @brief Valve (tanh) block.
 *
 * Models the pressure drop across a diode-like valve, which is implemented as a
 * non-linear piecewise resistor. See \cite Regazzoni2022
 * (equations 16 and 22).
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [D, l=$R_v$, *-*] (3,0)
 * node[anchor=south]{$P_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_{in}-P_{out}-Q_{in}\left[R(P_{out},P_{in})\right]=0
 * \f]
 *
 * \f[
 * Q_{in}-Q_{out}=0
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * R_i(p_1, p_2) =
 * \begin{cases}
 * R_{\min}, & p_1 < p_2, \\[0.5em]
 * R_{\max}, & p_1 \ge p_2.
 * \end{cases}
 * \f]
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{llll}P_{in} & Q_{in} &
 * P_{out} & Q_{out}\end{array}\right]^{T} \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{cccc}
 * 0 & 0 & 0 & 0 \\
 * 0 & 0 & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{cccc}
 * 1 & -(R(P_{in},P_{out})) & -1 & 0 \\
 * 0 &      1                 &  0 & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * 0 \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Rmax: Maximum (closed) valve resistance
 * * `1` Rmin: Minimum (open) valve resistance
 * * `2` upstream_block: Name of block connected upstream
 * * `3` downstream_block: Name of block connected downstream
 *
 */
class PiecewiseValve : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    RMAX = 0,
    RMIN = 1,
    STEEPNESS = 2,
  };

  /**
   * @brief Construct a new PiecewiseValve object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  PiecewiseValve(int id, Model* model)
      : Block(id, model, BlockType::piecewise_valve, BlockClass::valve,
              {{"Rmax", InputParameter()},
               {"Rmin", InputParameter()},
               {"upstream_block", InputParameter(false, false, false)},
               {"downstream_block", InputParameter(false, false, false)}}) {}

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set global_var_ids and global_eqn_ids of the element based on the
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
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{5, 0, 3};
};

#endif  // SVZERODSOLVER_MODEL_PiecewiseValve_HPP_
