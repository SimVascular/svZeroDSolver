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
 * @file junction.hpp
 * @brief MODEL::Junction source file
 */
#ifndef SVZERODSOLVER_MODEL_JUNCTION_HPP_
#define SVZERODSOLVER_MODEL_JUNCTION_HPP_

#include "../algebra/sparsesystem.hpp"
#include "block.hpp"

namespace MODEL {
/**
 * @brief Junction
 *
 * Models a junction with arbitrary inlets and outlets. Across all inlets and
 * outlets of the junction, mass is conserved and pressure is continuous.
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
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class Junction : public Block<T> {
 public:
  // Inherit constructors
  using Block<T>::Block;

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
  void update_constant(ALGEBRA::SparseSystem<T> &system,
                       std::vector<T> &parameters);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  std::map<std::string, int> num_triplets = {
      {"F", 0},
      {"E", 0},
      {"D", 0},
  };

  /**
   * @brief Get number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  std::map<std::string, int> get_num_triplets();

  unsigned int num_inlets;
  unsigned int num_outlets;
};

template <typename T>
void Junction<T>::setup_dofs(DOFHandler &dofhandler) {
  // Set number of equations of a junction block based on number of
  // inlets/outlets. Must be set before calling parent constructor
  num_inlets = this->inlet_nodes.size();
  num_outlets = this->outlet_nodes.size();
  Block<T>::setup_dofs_(dofhandler, num_inlets + num_outlets, {});
  num_triplets["F"] =
      (num_inlets + num_outlets - 1) * 2 + num_inlets + num_outlets;
}

template <typename T>
void Junction<T>::update_constant(ALGEBRA::SparseSystem<T> &system,
                                  std::vector<T> &parameters) {
  // Pressure conservation
  for (size_t i = 0; i < (num_inlets + num_outlets - 1); i++) {
    system.F.coeffRef(this->global_eqn_ids[i], this->global_var_ids[0]) = 1.0;
    system.F.coeffRef(this->global_eqn_ids[i],
                      this->global_var_ids[2 * i + 2]) = -1.0;
  }

  // Mass conservation
  for (size_t i = 1; i < num_inlets * 2; i = i + 2) {
    system.F.coeffRef(this->global_eqn_ids[num_inlets + num_outlets - 1],
                      this->global_var_ids[i]) = 1.0;
  }
  for (size_t i = (num_inlets * 2) + 1; i < (num_inlets + num_outlets) * 2;
       i = i + 2) {
    system.F.coeffRef(this->global_eqn_ids[num_inlets + num_outlets - 1],
                      this->global_var_ids[i]) = -1.0;
  }
}

template <typename T>
std::map<std::string, int> Junction<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_JUNCTION_HPP_