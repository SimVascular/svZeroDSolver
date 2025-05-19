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
 * @file ParallelRC.h
 * @brief model::ParallelRC source file
 */
#ifndef SVZERODSOLVER_MODEL_PARALLEL_RC_HPP_
#define SVZERODSOLVER_MODEL_PARALLEL_RC_HPP_

#include <math.h>

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Parallel RC circuit element for 0D blood flow modeling
 *
 * Models a parallel resistor and capacitor with shared pressure at
 * both the inlet and outlet.
 *
 * \f[
 * \begin{circuitikz} \draw
 *   node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 *   (1,0) node[anchor=south]{$P_{in}$}
 *   to [short, -*] (2,0)
 *   to [R, l=$R$] (5,0)
 *   to [short, -*] (6,0)
 *   node[anchor=south]{$P_{out}$}
 *   (7.2,0) -- (8,0) node[right] {$Q_{out}$};
 *   (2,0) to (2,-2) to [C, l=$C$] (6,-2) to (6,0);
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_\text{in}-P_\text{out} = 0 \f]
 *
 * \f[
 * Q_\text{in} - Q_R - Q_C - Q_\text{out} = 0 \f]
 * 
 * \f[
 * P_\text{in} - P_\text{out} - R \cdot Q_R = 0 \f]
 * 
 * \f[
 * Q_C - C (\dot{P}_\text{in} - \dot{P}_\text{out}) = 0 \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{llllll}P_{in} & Q_{in} &
 * P_{out} & Q_{out} & Q_R & Q_C\end{array}\right]^\text{T} \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Resistance value
 * * `1` Capacitance value
 *
 * ### Internal variables
 *
 * This block has two internal variables: Q_R and Q_C.
 *
 */
class ParallelRC : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    RESISTANCE = 0,
    CAPACITANCE = 1,
  };

  /**
   * @brief Construct a new ParallelRC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ParallelRC(int id, Model *model)
      : Block(id, model, BlockType::parallel_rc, BlockClass::vessel,
              {{"R", InputParameter()},
               {"C", InputParameter()}}) {
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
  TripletsContributions num_triplets{10, 1, 0};
};

#endif  // SVZERODSOLVER_MODEL_PARALLEL_RC_HPP_