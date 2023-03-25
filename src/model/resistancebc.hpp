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
 * @file resistancebc.hpp
 * @brief MODEL::ResistanceBC source file
 */
#ifndef SVZERODSOLVER_MODEL_RESISTANCEBC_HPP_
#define SVZERODSOLVER_MODEL_RESISTANCEBC_HPP_

#include "../algebra/sparsesystem.hpp"
#include "block.hpp"
#include "parameter.hpp"

namespace MODEL {

/**
 * @brief Resistance boundary condition.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1.0,0) to [R, l=$R$, *-*] (3,0)
 * node[anchor=south]{$P_{d}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P-P_d=R \cdot Q
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{ll}P^{e} & Q^{e}\end{array}\right]^{T}
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ll}1 & -R\end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{C}^{e}=\left[-P_d\right]
 * \f]
 *
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class ResistanceBC : public Block<T> {
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
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(ALGEBRA::SparseSystem<T> &system,
                   std::vector<T> &parameters);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  std::map<std::string, int> num_triplets = {
      {"F", 1},
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
};

template <typename T>
void ResistanceBC<T>::setup_dofs(DOFHandler &dofhandler) {
  Block<T>::setup_dofs_(dofhandler, 1, {});
}

template <typename T>
void ResistanceBC<T>::update_constant(ALGEBRA::SparseSystem<T> &system,
                                      std::vector<T> &parameters) {
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
}

template <typename T>
void ResistanceBC<T>::update_time(ALGEBRA::SparseSystem<T> &system,
                                  std::vector<T> &parameters) {
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
      -parameters[this->global_param_ids[0]];
  system.C(this->global_eqn_ids[0]) = -parameters[this->global_param_ids[1]];
}

template <typename T>
std::map<std::string, int> ResistanceBC<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_RESISTANCEBC_HPP_
