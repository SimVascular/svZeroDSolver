// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
 * non-linear hyperbolic-tangent resistor.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1.0,0) to [D, l=$R_v$, *-*] (3,0)
 * node[anchor=south]{$P_{d}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_{in}-P_{out}-Q_{in}\left[R_{min} +
 * (R_{max}-R_{min})\frac{1}{2}(1+tanh\{k(P_{out}-P_P{in})\})\right]=0 \f]
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
 * \frac{1}{2} k Q_{in} (R_{max}-R_{min})
 * \left[1-tanh^2\{k(P_{out}-P_{in})\}\right] &
 * -\frac{1}{2}(R_{max}-R_{min})tanh\{k(P_{out}-P_{in})\} &
 * -\frac{1}{2} k Q_{in} (R_{max}-R_{min})
 * \left[1-tanh^2\{k(P_{out}-P_{in})\}\right] & 0
 * \\ 0 & 0 & 0 & 0 \end{array}\right] \f]
 *
 * \f[
 * \left(\frac{\partial\mathbf{c}}{\partial\dot{\mathbf{y}}}\right)^{e} =
 * \left[\begin{array}{cccc}
 * 0 & 0 & 0 & 0 \\
 * 0 & 0 & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * See \cite pfaller2019importance.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Maximum (closed) valve resistance
 * * `1` Minimum (open) valve resistance
 * * `2` Steepness of sigmoid function
 *
 */
class ValveTanh : public Block {
 public:
  // Inherit constructors
  using Block::Block;

  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    RMAX = 0,
    RMIN = 1,
    STEEPNESS = 2,
  };

  // explicit ValveTanh(int id, const std::vector<int> &param_ids, Model *model)
  //    : Block(id, param_ids, model){};

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
