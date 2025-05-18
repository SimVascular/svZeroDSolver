// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ValveTanh.h
 * @brief model::ValveTanh source file
 */
#ifndef SVZERODSOLVER_MODEL_VALVETANH_HPP_
#define SVZERODSOLVER_MODEL_VALVETANH_HPP_

#include <math.h>

#include "Block.h"
#include "SparseSystem.h"
#include "debug.h"

/**
 * @brief Valve (tanh) block.
 *
 * Models the pressure drop across a diode-like valve, which is implemented as a
 * non-linear hyperbolic-tangent resistor. See \cite pfaller2019importance
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
 * P_{in}-P_{out}-Q_{in}\left[R_{min} +
 * (R_{max}-R_{min})\frac{1}{2}\left[1+tanh\{k(P_{out}-P{in})\}\right]\right]=0
 * \f]
 *
 * \f[
 * Q_{in}-Q_{out}=0
 * \f]
 *
 * ### Local contributions
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
 * 1 & -(R_{max}+R_{min})/2.0 & -1 & 0 \\
 * 0 &      1                 &  0 & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * -\frac{1}{2}Q_{in}(R_{max}-R_{min})tanh\{k(P_{out}-P_{in})\} \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \left(\frac{\partial\mathbf{c}}{\partial\mathbf{y}}\right)^{e} =
 * \left[\begin{array}{cccc}
 * A & B & C & 0 \\
 * 0 & 0 & 0 & 0 \end{array}\right] \f]
 * where,
 * \f[
 * A = \frac{1}{2} k Q_{in}
 * (R_{max}-R_{min})\left[1-tanh^2\{k(P_{out}-P_{in})\}\right] \\
 * \f]
 * \f[
 * B = -\frac{1}{2}(R_{max}-R_{min})tanh\{k(P_{out}-P_{in})\} \\
 * \f]
 * \f[
 * C = -\frac{1}{2} k Q_{in}
 * (R_{max}-R_{min})\left[1-tanh^2\{k(P_{out}-P_{in})\}\right] \f]
 *
 * \f[
 * \left(\frac{\partial\mathbf{c}}{\partial\dot{\mathbf{y}}}\right)^{e} =
 * \left[\begin{array}{cccc}
 * 0 & 0 & 0 & 0 \\
 * 0 & 0 & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Rmax: Maximum (closed) valve resistance
 * * `1` Rmin: Minimum (open) valve resistance
 * * `2` Steepness: Steepness of sigmoid function
 * * `3` upstream_block: Name of block connected upstream
 * * `4` downstream_block: Name of block connected downstream
 *
 * ### Internal variables
 *
 * This block has no internal variables.
 *
 */
class ValveTanh : public Block {
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
   * @brief Construct a new ValveTanh object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ValveTanh(int id, Model *model)
      : Block(id, model, BlockType::valve_tanh, BlockClass::valve,
              {{"Rmax", InputParameter()},
               {"Rmin", InputParameter()},
               {"Steepness", InputParameter()},
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
  TripletsContributions num_triplets{5, 0, 3};
};

#endif  // SVZERODSOLVER_MODEL_VALVETANH_HPP_
