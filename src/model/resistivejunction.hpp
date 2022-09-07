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
 * @file resistivejunction.hpp
 * @brief MODEL::ResistiveJunction source file
 */
#ifndef SVZERODSOLVER_MODEL_RESISTIVEJUNCTION_HPP_
#define SVZERODSOLVER_MODEL_RESISTIVEJUNCTION_HPP_

#include "../algebra/sparsesystem.hpp"
#include "block.hpp"

namespace MODEL {
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
 * (3,0) node[anchor=south]{} to [R, l=$R_{out,1}$, -*] (5,1.0); \draw (4.3,1.1)
 * node[anchor=south] {$P_{out,1}$}; \draw (3,0) node[anchor=south]{} to [R,
 * l=$R_{out,2}$, -*] (5,-1.0); \draw (4.3,-1.1) node[anchor=north]
 * {$P_{out,2}$}; \draw [-latex] (5.15,1.1) -- (5.75,1.4) node[right]
 * {$Q_{out,1}$}; \draw [-latex] (5.15,-1.1) -- (5.75,-1.4) node[right]
 * {$Q_{out,2}$}; \end{circuitikz} \f]
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
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class ResistiveJunction : public Block<T> {
 public:
  /**
   * @brief Parameters of the element.
   *
   * Struct containing all constant and/or time-dependent parameters of the
   * element.
   */
  struct Parameters : public Block<T>::Parameters {
    std::vector<T> R;  ///< Poiseuille resistance
  };

  /**
   * @brief Construct a new ResistiveJunction object
   *
   * @param name Name
   */
  ResistiveJunction(std::vector<T> R, std::string name);

  /**
   * @brief Destroy the ResistiveJunction object
   *
   */
  ~ResistiveJunction();

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
   */
  void update_constant(ALGEBRA::SparseSystem<T> &system);

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

 private:
  Parameters params;
  unsigned int num_inlets;
  unsigned int num_outlets;
};

template <typename T>
ResistiveJunction<T>::ResistiveJunction(std::vector<T> R, std::string name)
    : Block<T>(name) {
  this->name = name;
  this->params.R = R;
}

template <typename T>
ResistiveJunction<T>::~ResistiveJunction() {}

template <typename T>
void ResistiveJunction<T>::setup_dofs(DOFHandler &dofhandler) {
  // Set number of equations of a junction block based on number of
  // inlets/outlets. Must be set before calling parent constructor
  num_inlets = this->inlet_nodes.size();
  num_outlets = this->outlet_nodes.size();
  Block<T>::setup_dofs_(dofhandler, num_inlets + num_outlets + 1,
                        {"pressure_c"});
  num_triplets["F"] = (num_inlets + num_outlets) * 4;
}

template <typename T>
void ResistiveJunction<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  for (size_t i = 0; i < num_inlets; i++) {
    system.F.coeffRef(this->global_eqn_ids[i], this->global_var_ids[i * 2]) =
        1.0;
    system.F.coeffRef(this->global_eqn_ids[i],
                      this->global_var_ids[i * 2 + 1]) = -params.R[i];
    system.F.coeffRef(this->global_eqn_ids[i], this->global_var_ids.back()) =
        -1.0;
  }
  for (size_t i = num_inlets; i < num_inlets + num_outlets; i++) {
    system.F.coeffRef(this->global_eqn_ids[i], this->global_var_ids[i * 2]) =
        -1.0;
    system.F.coeffRef(this->global_eqn_ids[i],
                      this->global_var_ids[i * 2 + 1]) = -params.R[i];
    system.F.coeffRef(this->global_eqn_ids[i], this->global_var_ids.back()) =
        1.0;
  }

  // Mass conservation
  for (size_t i = 1; i < num_inlets * 2; i = i + 2) {
    system.F.coeffRef(this->global_eqn_ids[num_inlets + num_outlets],
                      this->global_var_ids[i]) = 1.0;
  }
  for (size_t i = (num_inlets * 2) + 1; i < (num_inlets + num_outlets) * 2;
       i = i + 2) {
    system.F.coeffRef(this->global_eqn_ids[num_inlets + num_outlets],
                      this->global_var_ids[i]) = -1.0;
  }
}

template <typename T>
std::map<std::string, int> ResistiveJunction<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_RESISTIVEJUNCTION_HPP_