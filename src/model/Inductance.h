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
 * @file Inductance.h
 * @brief model::Inductance source file
 */
#ifndef SVZERODSOLVER_MODEL_INDUCTANCE_HPP_
#define SVZERODSOLVER_MODEL_INDUCTANCE_HPP_

#include <math.h>

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Simple inductance element for 0D blood flow modeling
 *
 * Models a pure inductance element relating pressure drop to change in flow.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [L, l=$L$, *-*] (7,0)
 * node[anchor=south]{$P_{out}$};
 * \draw [-latex] (7.2,0) -- (8,0) node[right] {$Q_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_\text{in}-P_\text{out} - L \dot{Q}_\text{out} = 0 \f]
 *
 * \f[
 * Q_\text{in}-Q_\text{out} = 0 \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{llll}P_{i n} & Q_{in} &
 * P_{out} & Q_{out}\end{array}\right]^\text{T} \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{cccc}
 * 1 & 0 & -1 &  0 \\
 * 0 & 1 &  0 & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{cccc}
 *  0 &  0 & 0 & -L \\
 *  0 &  0 & 0 &  0
 * \end{array}\right]
 * \f]
 *
 * ### Gradient
 *
 * Gradient of the equations with respect to the parameters:
 *
 * \f[
 * \mathbf{J}^{e} = \left[\begin{array}{cccc}
 * -\dot{y}_4 \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Inductance value
 *
 * ### Internal variables
 *
 * This block has no internal variables.
 *
 */
class Inductance : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    INDUCTANCE = 0,
  };

  /**
   * @brief Construct a new Inductance object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  Inductance(int id, Model *model)
      : Block(id, model, BlockType::inductance, BlockClass::vessel,
              {{"L", InputParameter()}}) {
    is_vessel = true;
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
   * @brief Update the constant contributions of the element in a sparse
   * system
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
  TripletsContributions num_triplets{4, 1, 0};
};

#endif  // SVZERODSOLVER_MODEL_INDUCTANCE_HPP_