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
#include "blocktype.hpp"
#include "bloodvessel.hpp"

namespace MODEL {
/**
 * @brief BloodVesselJunction
 *
 * Models a junction with one inlet and arbitrary outlets using
 * modified blood vessel elements between each inlet and outlet pair.
 *
 * \f[
 * \begin{circuitikz}
 * \draw node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0.1) node[anchor=south]{$P_{in}$};
 * \draw (1,0) to [short, *-] (2.5,0.75);
 * \draw (1,0) to [short, *-] (2.5,-0.75);
 * \draw (2.5,0.75) node[anchor=south]{} to [generic, l_=$BV_{1}$, -*]
 * (4.5,0.75); \draw (2.4,0.75) node[anchor=south]{}; \draw (4.6,0.75)
 * node[anchor=south] {$P_{out,1}$}; \draw (2.5,-0.75) node[anchor=south]{} to
 * [generic, l^=$BV_{2}$, -*] (4.5,-0.75); \draw (2.4,-0.75)
 * node[anchor=north]{}; \draw (4.6,-0.75) node[anchor=north]
 * {$P_{out,2}$}; \draw [-latex] (4.7,0.75) -- (5.5,0.75) node[right]
 * {$Q_{out,1}$}; \draw [-latex] (4.7,-0.75) -- (5.5,-0.75) node[right]
 * {$Q_{out,2}$}; \end{circuitikz} \f]
 *
 * Each blood vessel is modelled as:
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R$, *-] (3,0)
 * to [R, l=$R_{ste}$, -] (5,0)
 * (5,0) to [L, l=$L$, -*] (7,0)
 * node[anchor=south]{$P_{out}$};
 * \draw [-latex] (7.2,0) -- (8,0) node[right] {$Q_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * Q_{in}-\sum_{i}^{n_{outlets}} Q_{out, i}
 * \f]
 *
 * \f[
 * P_{in}-P_{out,i} - (R+R_{ste}) \cdot Q_{out,i} -
 * L \frac{d Q_{out,i}}{d t} \quad \forall i \in n_{outlets} \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lllllll}P_{in}^{e} & Q_{in}^{e}
 * & P_{out, 1}^{e} & Q_{out, 1}^{e} &
 * \dots & P_{out, i}^{e} & Q_{out, i}^{e}\end{array}\right] \f]
 *
 * \f[
 * \mathbf{F}^{e} = \left[\begin{array}{lllllllll}
 * 0 & 1 & 0 & -1 & 0 & -1 & 0 & -1 & \dots \\
 * 1 & 0 & -1 & -R_{1}-R_{ste,1} & 0 & 0 & 0 & 0 & \dots\\
 * 1 & 0 & 0 & 0 & -1 & -R_{2}-R_{ste,2} & 0& 0& \dots \\
 * \vdots & & & & & \ddots & \ddots & &
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{E}^{e} = \left[\begin{array}{lllllllll}
 * 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \dots \\
 * 0 & 0 & 0 & -L_{1} & 0 & 0 & 0 & 0 & \dots\\
 * 0 & 0 & 0 & 0 & 0 & -L_{2} & 0 & 0 & \dots\\
 * & & & & & \ddots & \ddots & &
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{D}^{e} = \left[\begin{array}{lllllllll}
 * 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \dots \\
 * 0 & 0 & 0 &
 * -R_{ste,1} & 0 & 0 & 0 & 0 & \dots\\
 * 0 & 0 & 0 & 0 & 0 & -R_{ste,2} & 0 & 0 & \dots\\ & & & & & \ddots & \ddots &
 * & \end{array}\right] \f]
 *
 * ### Gradient
 *
 * Gradient of the equations with respect to the parameters:
 *
 * \f[
 * \mathbf{J}^{e} = \left[\begin{array}{lllllllll}
 * 0 & 0 & \dots & 0 & 0 & \dots & 0 & 0 &
 * \dots \\
 * - y_4 & 0 & \dots & - \dot y_4 & 0 & \dots & |y_4| y_4 & 0 & \dots \\
 * 0 & - y_6 & \dots & 0 & - \dot y_6 & \dots & 0 & |y_6| y_6 & \dots \\
 * 0 & 0 & \ddots & 0 & 0 & \ddots & 0 & 0 & \ddots \\
 * \end{array}\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `i` Poiseuille resistance for inner blood vessel `i`
 * * `i+num_outlets` Inductance for inner blood vessel `i`
 * * `i+2*num_outlets` Stenosis coefficient for inner blood vessel `i`
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class BloodVesselJunction : public Block<T> {
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
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  virtual void update_solution(ALGEBRA::SparseSystem<T> &system,
                               std::vector<T> &parameters,
                               Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
                               Eigen::Matrix<T, Eigen::Dynamic, 1> &dy);

  /**
   * @brief Set the gradient of the block contributions with respect to the
   * parameters
   *
   * @param jacobian Jacobian with respect to the parameters
   * @param alpha Current parameter vector
   * @param residual Residual with respect to the parameters
   * @param y Current solution
   * @param dy Time-derivative of the current solution
   */
  void update_gradient(Eigen::SparseMatrix<T> &jacobian,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &residual,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &alpha,
                       std::vector<T> &y, std::vector<T> &dy);

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
  int num_outlets;
};

template <typename T>
void BloodVesselJunction<T>::setup_dofs(DOFHandler &dofhandler) {
  if (this->inlet_nodes.size() != 1) {
    throw std::runtime_error(
        "Blood vessel junction does not support multiple inlets.");
  }
  num_outlets = this->outlet_nodes.size();
  Block<T>::setup_dofs_(dofhandler, num_outlets + 1, {});
  num_triplets["F"] = 1 + 4 * num_outlets;
  num_triplets["E"] = 3 * num_outlets;
  num_triplets["D"] = 2 * num_outlets;
}

template <typename T>
void BloodVesselJunction<T>::update_constant(ALGEBRA::SparseSystem<T> &system,
                                             std::vector<T> &parameters) {
  // Mass conservation
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = 1.0;
  for (size_t i = 0; i < num_outlets; i++) {
    T inductance = parameters[this->global_param_ids[num_outlets + i]];
    system.F.coeffRef(this->global_eqn_ids[0],
                      this->global_var_ids[3 + 2 * i]) = -1.0;

    system.F.coeffRef(this->global_eqn_ids[i + 1], this->global_var_ids[0]) =
        1.0;
    system.F.coeffRef(this->global_eqn_ids[i + 1],
                      this->global_var_ids[2 + 2 * i]) = -1.0;

    system.E.coeffRef(this->global_eqn_ids[i + 1],
                      this->global_var_ids[3 + 2 * i]) = -inductance;
  }
}

template <typename T>
void BloodVesselJunction<T>::update_solution(
    ALGEBRA::SparseSystem<T> &system, std::vector<T> &parameters,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &dy) {
  for (size_t i = 0; i < num_outlets; i++) {
    // Get parameters
    T resistance = parameters[this->global_param_ids[i]];
    T stenosis_coeff = parameters[this->global_param_ids[2 * num_outlets + i]];
    T q_out = y[this->global_var_ids[3 + 2 * i]];
    T stenosis_resistance = stenosis_coeff * fabs(q_out);

    // Mass conservation
    system.F.coeffRef(this->global_eqn_ids[i + 1],
                      this->global_var_ids[3 + 2 * i]) =
        -resistance - stenosis_resistance;

    system.D.coeffRef(this->global_eqn_ids[i + 1],
                      this->global_var_ids[3 + 2 * i]) = -stenosis_resistance;
  }
}

template <typename T>
void BloodVesselJunction<T>::update_gradient(
    Eigen::SparseMatrix<T> &jacobian,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &residual,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &alpha, std::vector<T> &y,
    std::vector<T> &dy) {
  T p_in = y[this->global_var_ids[0]];
  T q_in = y[this->global_var_ids[1]];

  residual(this->global_eqn_ids[0]) = q_in;
  for (size_t i = 0; i < num_outlets; i++) {
    // Get parameters
    T resistance = alpha[this->global_param_ids[i]];
    T inductance = alpha[this->global_param_ids[num_outlets + i]];
    T stenosis_coeff = 0.0;
    if (this->global_param_ids.size() / num_outlets > 2) {
      stenosis_coeff = alpha[this->global_param_ids[2 * num_outlets + i]];
    }
    T q_out = y[this->global_var_ids[3 + 2 * i]];
    T p_out = y[this->global_var_ids[2 + 2 * i]];
    T dq_out = dy[this->global_var_ids[3 + 2 * i]];
    T stenosis_resistance = stenosis_coeff * fabs(q_out);

    // Resistance
    jacobian.coeffRef(this->global_eqn_ids[i + 1], this->global_param_ids[i]) =
        -q_out;

    // Inductance
    jacobian.coeffRef(this->global_eqn_ids[i + 1],
                      this->global_param_ids[num_outlets + i]) = -dq_out;

    // Stenosis Coefficent
    if (this->global_param_ids.size() / num_outlets > 2) {
      jacobian.coeffRef(this->global_eqn_ids[i + 1],
                        this->global_param_ids[2 * num_outlets + i]) =
          -fabs(q_out) * q_out;
    }

    residual(this->global_eqn_ids[0]) -= q_out;
    residual(this->global_eqn_ids[i + 1]) =
        p_in - p_out - (resistance + stenosis_resistance) * q_out -
        inductance * dq_out;
  }
}

template <typename T>
std::map<std::string, int> BloodVesselJunction<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_BLOODVESSELJUNCTION_HPP_