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
  /**
   * @brief Parameters of the element.
   *
   * Struct containing all constant and/or time-dependent parameters of the
   * element.
   */
  struct Parameters : public Block<T>::Parameters {
    Parameter<T> R;   ///< Time-dependent resistance
    Parameter<T> Pd;  ///< Time-dependent distal pressure
  };

  /**
   * @brief Construct a new ResistanceBC object
   *
   * @param R Time-dependent resistance
   * @param Pd Time-dependent distal pressure
   * @param name Name
   */
  ResistanceBC(Parameter<T> R, Parameter<T> Pd, std::string name);

  /**
   * @brief Destroy the ResistanceBC object
   *
   */
  ~ResistanceBC();

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
   * @brief Update parameters of a block.
   *
   * @param params New parameters.
   */
  void update_block_params(std::vector<T> new_params);

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   */
  void update_constant(ALGEBRA::SparseSystem<T> &system);

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param time Current time
   */
  void update_time(ALGEBRA::SparseSystem<T> &system, T time);

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

  /**
   * @brief Convert the block to a steady behavior
   *
   * Converts the resistance and distal pressure to the constant means of
   * themselve
   */
  void to_steady();

  /**
   * @brief Convert the block to a steady behavior
   *
   */
  void to_unsteady();

 private:
  Parameters params;
};

template <typename T>
ResistanceBC<T>::ResistanceBC(Parameter<T> R, Parameter<T> Pd, std::string name)
    : Block<T>(name) {
  this->name = name;
  this->params.R = R;
  this->params.Pd = Pd;
}

template <typename T>
ResistanceBC<T>::~ResistanceBC() {}

template <typename T>
void ResistanceBC<T>::setup_dofs(DOFHandler &dofhandler) {
  Block<T>::setup_dofs_(dofhandler, 1, {});
}

template <typename T>
void ResistanceBC<T>::update_block_params(std::vector<T> new_params) {
  std::vector<T> t_new;
  std::vector<T> R_new;
  std::vector<T> Pd_new;
  int num_time_pts = (int)new_params[0];
  for (int i = 0; i < num_time_pts; i++) {
    t_new.push_back(new_params[1 + i]);
    R_new.push_back(new_params[1 + num_time_pts + i]);
    Pd_new.push_back(new_params[1 + 2 * num_time_pts + i]);
  }
  this->params.R.update(t_new, R_new);
  this->params.Pd.update(t_new, Pd_new);
}
template <typename T>
void ResistanceBC<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
}

template <typename T>
void ResistanceBC<T>::update_time(ALGEBRA::SparseSystem<T> &system, T time) {
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
      -params.R.get(time);
  system.C(this->global_eqn_ids[0]) = -params.Pd.get(time);
}

template <typename T>
void ResistanceBC<T>::to_steady() {
  params.R.to_steady();
  params.Pd.to_steady();
}

template <typename T>
void ResistanceBC<T>::to_unsteady() {
  params.R.to_unsteady();
  params.Pd.to_unsteady();
}

template <typename T>
std::map<std::string, int> ResistanceBC<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_RESISTANCEBC_HPP_
