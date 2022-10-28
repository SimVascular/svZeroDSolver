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
 * @file flowreferencebc.hpp
 * @brief MODEL::FlowReferenceBC source file
 */
#ifndef SVZERODSOLVER_MODEL_FLOWREFERENCEBC_HPP_
#define SVZERODSOLVER_MODEL_FLOWREFERENCEBC_HPP_

#include "../algebra/sparsesystem.hpp"
#include "block.hpp"
#include "timedependentparameter.hpp"

namespace MODEL {

/**
 * @brief Flow reference boundary condition.
 *
 * Applies a prescribed flow to a boundary.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$\hat{Q}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P$} to [short, *-] (1.2,0) ;
 * \draw [-latex] (1.4,0) -- (2.2,0) node[right] {$Q$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * Q=\hat{Q}
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{ll}P^{e} & Q^{e}\end{array}\right]^{T}
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ll}0 & 1\end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{C}^{e}=\left[\hat{Q}\right]
 * \f]
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class FlowReferenceBC : public Block<T> {
 public:
  /**
   * @brief Parameters of the element.
   *
   * Struct containing all constant and/or time-dependent parameters of the
   * element.
   */
  struct Parameters : public Block<T>::Parameters {
    TimeDependentParameter<T> Q;  ///< Time-dependent flow
  };

  /**
   * @brief Construct a new FlowReferenceBC object
   *
   * @param Q Time dependent flow
   * @param name Name
   */
  //FlowReferenceBC(TimeDependentParameter<T> Q, std::string name);
  FlowReferenceBC(TimeDependentParameter<T> Q, std::string name, std::string coupling_loc = "None");

  /**
   * @brief Destroy the FlowReferenceBC object
   *
   */
  ~FlowReferenceBC();

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
   * Converts the prescribed flow to the constant mean of itself
   *
   */
  void to_steady();

  /**
   * @brief Convert the block to an unsteady behavior
   *
   */
  void to_unsteady();

  /**
   * @brief Specify is this is an inlet or outlet to the svZeroD model when used for external coupling.
   *
   */
  std::string coupling_loc;

 private:
  Parameters params;
  bool external_coupling = false;
};

template <typename T>
FlowReferenceBC<T>::FlowReferenceBC(TimeDependentParameter<T> Q,
                                    std::string name, std::string coupling_loc)
    : Block<T>(name) {
  this->name = name;
  this->params.Q = Q;
  this->coupling_loc = coupling_loc;
  if (coupling_loc != "None") {
    this->external_coupling = true;
  }
}

template <typename T>
FlowReferenceBC<T>::~FlowReferenceBC() {}

template <typename T>
void FlowReferenceBC<T>::setup_dofs(DOFHandler &dofhandler) {
  Block<T>::setup_dofs_(dofhandler, 1, {});
}

template <typename T>
void FlowReferenceBC<T>::update_block_params(std::vector<T> new_params) {
  std::vector<T> t_new;
  std::vector<T> Q_new;
  //std::cout << "[update_block_params] new_params: ";
  int num_time_pts = (int) new_params[0];
  for (int i = 0; i < num_time_pts; i++) {
    t_new.push_back(new_params[1+i]);
    Q_new.push_back(new_params[1+num_time_pts+i]);
  }
  this->params.Q.update_params(t_new,Q_new);
  std::cout << "[update_block_params] params.Q.values: " << std::endl;
  for(int j = 0; j < 2; j++) {
    std::cout << params.Q.values[j] << " ";
  }
  std::cout<<std::endl;
}

template <typename T>
void FlowReferenceBC<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = 1.0;
}

template <typename T>
void FlowReferenceBC<T>::update_time(ALGEBRA::SparseSystem<T> &system, T time) {
  system.C(this->global_eqn_ids[0]) = -params.Q.get(time);
//std::cout<<"[FlowReferenceBC<T>::update_time] name, time: "<<this->name<<", "<<time<<std::endl;
//std::cout<<"params.Q.get(time) = "<<params.Q.get(time, this->external_coupling)<<std::endl;
}

template <typename T>
void FlowReferenceBC<T>::to_steady() {
  params.Q.to_steady();
}

template <typename T>
void FlowReferenceBC<T>::to_unsteady() {
  params.Q.to_unsteady();
}

template <typename T>
std::map<std::string, int> FlowReferenceBC<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_FLOWREFERENCEBC_HPP_
