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
 * @file WindkesselBC.h
 * @brief model::WindkesselBC source file
 */
#ifndef SVZERODSOLVER_MODEL_WINDKESSELBC_HPP_
#define SVZERODSOLVER_MODEL_WINDKESSELBC_HPP_

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Windkessel RCR boundary condition.
 *
 * Models the mechanical behavior of a Windkessel boundary condition.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_p$, *-] (3,0)
 * node[anchor=south]{$P_{c}$}
 * to [R, l=$R_d$, *-*] (5,0)
 * node[anchor=south]{$P_{ref}$}
 * (3,0) to [C, l=$C$, *-] (3,-1.5)
 * node[ground]{};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * R_{d} Q^{e}-P_{c}^{e}+P_{r e f}-R_{d} C \frac{d P_{c}^{e}}{d t}=0
 * \f]
 *
 * \f[
 * P^{e}-P_{c}^{e}-R_{p} Q^{e}=0
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lll}P^{e} & Q^{e} &
 * P_{c}^{e}\end{array}\right]^{T} \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{ccc}
 * 0 & 0 & -R_{d} C \\
 * 0 & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccc}
 * 0 & R_{d} & -1 \\
 * 1 & -R_{p} & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * P_{r e f} \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * See \cite vignon04.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Proximal resistance
 * * `1` Capacitance
 * * `2` Distal resistance
 * * `3` Distal pressure
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `pressure_c`: Pressure at the capacitor
 *
 */
class WindkesselBC : public Block {
 public:
  /**
   * @brief Construct a new WindkesselBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  WindkesselBC(int id, Model *model)
      : Block(id, model, BlockType::windkessel_bc,
              BlockClass::boundary_condition,
              {{"Rp", InputParameter()},
               {"C", InputParameter()},
               {"Rd", InputParameter()},
               {"Pd", InputParameter(true)}}) {}

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
  TripletsContributions num_triplets{5, 1, 0};
};

#endif  // SVZERODSOLVER_MODEL_WINDKESSELBC_HPP_
