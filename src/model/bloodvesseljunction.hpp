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
 * @file bloodvesseljunction.hpp
 * @brief MODEL::BloodVesselJunction source file
 */
#ifndef SVZERODSOLVER_MODEL_BLOODVESSELJUNCTION_HPP_
#define SVZERODSOLVER_MODEL_BLOODVESSELJUNCTION_HPP_

#include "../algebra/sparsesystem.hpp"
#include "block.hpp"
#include "bloodvessel.hpp"

namespace MODEL {
/**
 * @brief BloodVesselJunction
 *
 * Models a junction with one inlet and arbitrary outlets using
 * MODEL::BloodVessel elements between each inlet and outlet pair.
 *
 * \f[
 * \begin{circuitikz}
 * \draw node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0.1) node[anchor=south]{$P_{in}$};
 * \draw (1,0) to [short, *-] (2.5,0.75);
 * \draw (1,0) to [short, *-] (2.5,-0.75);
 * \draw (2.5,0.75) node[anchor=south]{} to [generic, l_=$BV_{1}$, -*]
 * (4.5,0.75); \draw (2.4,0.75) node[anchor=south]{$Q_{in,1}$}; \draw (4.6,0.75)
 * node[anchor=south] {$P_{out,1}$}; \draw (2.5,-0.75) node[anchor=south]{} to
 * [generic, l^=$BV_{2}$, -*] (4.5,-0.75); \draw (2.4,-0.75)
 * node[anchor=north]{$Q_{in,2}$}; \draw (4.6,-0.75) node[anchor=north]
 * {$P_{out,2}$}; \draw [-latex] (4.7,0.75) -- (5.5,0.75) node[right]
 * {$Q_{out,1}$}; \draw [-latex] (4.7,-0.75) -- (5.5,-0.75) node[right]
 * {$Q_{out,2}$}; \end{circuitikz} \f]
 *
 * ### Governing equations
 *
 * The governing equations are mainly defined by the blood vessel elements
 * (see MODEL::BloodVessel for more information). For the sake of brevity,
 * the blood vessel contributions are not repeated here. In addition to the
 * individual blood vessel contributions, the following equation ensures
 * mass conservation between the inlets of the blood vessel elements:
 *
 * \f[
 * Q_{in}=\sum_{i}^{n_{outlets}} Q_{in, i}
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{llllllllll}P_{in, 1}^{e} & Q_{in, 1}^{e}
 * & P_{out, 1}^{e} & Q_{out, 1}^{e} &
 * \dots & P_{out, i}^{e} & Q_{out, i}^{e} & Q_{in,1} & \dots &
 * Q_{in,i}\end{array}\right] \f]
 *
 * Mass conservation
 *
 * \f[
 * \mathbf{F}^{e} = \left[\begin{array}{llllllllll}
 * 0 & 1 & 0 & 0 & \dots & 0 & 0 & -1 & \dots & -1 \\
 * & & & & & & & & & \\
 * \multicolumn{10}{c}{\mathrm{< blood \, vessel \, contributions >}} \\
 * & & & & & & & & &
 * \end{array}\right]
 * \f]
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class BloodVesselJunction : public Block<T> {
 public:
  /**
   * @brief Parameters of the element.
   *
   * Struct containing all constant and/or time-dependent parameters of the
   * element.
   */
  struct Parameters : public Block<T>::Parameters {
    std::vector<T> R;                     ///< Poiseuille resistances
    std::vector<T> C;                     ///< Capacitances
    std::vector<T> L;                     ///< Inductances
    std::vector<T> stenosis_coefficient;  ///< Stenosis Coefficients
  };

  /**
   * @brief Construct a new BloodVesselJunction object
   *
   * @param R Poiseuille resistances
   * @param C Capacitances
   * @param L Inductances
   * @param stenosis_coefficient Stenosis Coefficients
   * @param name Name
   */
  BloodVesselJunction(std::vector<T> R, std::vector<T> C, std::vector<T> L,
                      std::vector<T> stenosis_coefficient, std::string name);

  /**
   * @brief Destroy the BloodVesselJunction object
   *
   */
  ~BloodVesselJunction();

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
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param y Current solution
   */
  void update_solution(ALGEBRA::SparseSystem<T> &system,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

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
  std::vector<BloodVessel<T> *> blood_vessels;
  Parameters params;
  unsigned int num_outlets;
};

template <typename T>
BloodVesselJunction<T>::BloodVesselJunction(std::vector<T> R, std::vector<T> C,
                                            std::vector<T> L,
                                            std::vector<T> stenosis_coefficient,
                                            std::string name)
    : Block<T>(name) {
  this->name = name;
  this->params.R = R;
  this->params.C = C;
  this->params.L = L;
  this->params.stenosis_coefficient = stenosis_coefficient;
}

template <typename T>
BloodVesselJunction<T>::~BloodVesselJunction() {}

template <typename T>
void BloodVesselJunction<T>::setup_dofs(DOFHandler &dofhandler) {
  // Set number of equations of a junction block based on number of
  // inlets/outlets. Must be set before calling parent constructor
  num_outlets = this->outlet_nodes.size();
  std::list<std::string> internal_var_names;
  for (size_t i = 0; i < num_outlets; i++) {
    internal_var_names.push_back("flow_" + std::to_string(i));
  }
  Block<T>::setup_dofs_(dofhandler, 1, internal_var_names);
  for (size_t i = 0; i < num_outlets; i++) {
    blood_vessels.push_back(new BloodVessel<T>(
        params.R[i], params.C[i], params.L[i], params.stenosis_coefficient[i],
        this->name + "_bv" + std::to_string(i)));
    blood_vessels[i]->inlet_nodes.push_back(this->inlet_nodes[0]);
    blood_vessels[i]->outlet_nodes.push_back(this->outlet_nodes[i]);
    blood_vessels[i]->setup_dofs(dofhandler);
    blood_vessels[i]->global_var_ids[1] =
        this->global_var_ids.end()[i - num_outlets];
  }
  num_triplets = {
      {"F", 10 * num_outlets + num_outlets + 1},
      {"E", 2 * num_outlets},
      {"D", 2 * num_outlets},
  };
}

template <typename T>
void BloodVesselJunction<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  for (auto bv : blood_vessels) {
    bv->update_constant(system);
  }
  // Mass conservation
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = 1.0;
  for (size_t i = 0; i < num_outlets; i++) {
    system.F.coeffRef(this->global_eqn_ids[0],
                      this->global_var_ids.end()[-(i + 1)]) = -1.0;
  }
}

template <typename T>
void BloodVesselJunction<T>::update_solution(
    ALGEBRA::SparseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  for (auto bv : blood_vessels) {
    bv->update_solution(system, y);
  }
}

template <typename T>
std::map<std::string, int> BloodVesselJunction<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_BLOODVESSELJUNCTION_HPP_