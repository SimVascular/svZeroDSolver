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
 * @file bloodvessel.hpp
 * @brief MODEL::BloodVessel source file
 */
#ifndef SVZERODSOLVER_MODEL_BLOODVESSEL_HPP_
#define SVZERODSOLVER_MODEL_BLOODVESSEL_HPP_

#include <math.h>

#include "../algebra/sparsesystem.hpp"
#include "block.hpp"

namespace MODEL {

/**
 * @brief Resistor-capacitor-inductor blood vessel with optional stenosis
 *
 * Models the mechanical behavior of a bloodvessel with optional stenosis.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R$, *-] (3,0)
 * to [R, l=$R_{ste}$, -] (5,0) node[anchor=south] {$P_{C}$}
 * (5,0) to [L, l=$L$, *-*] (7,0)
 * node[anchor=south]{$P_{out}$}
 * (5,0) to [C, l=$C$, *-] (5,-1.5)
 * node[ground]{};
 * \draw [-latex] (7.2,0) -- (8,0) node[right] {$Q_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_{in}^{e}-P_{out}^{e}-(R+R_{ste}) Q_{in}^{e}-L\frac{d Q_{out}^{e}}{dt}=0
 * \f]
 *
 * \f[
 * Q_{i n}^{e}-Q_{o u t}^{e}-C \frac{d P_{c}^{e}}{d t}=0
 * \f]
 *
 * \f[
 * P_{i n}^{e}-(R+R_{ste}) Q_{i n}^{e}-P_{c}=0
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lllll}P_{i n}^{e} & Q_{in}^{e} &
 * P_{out}^{e} & Q_{out}^{e} & P_C\end{array}\right]^{T} \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{ccccc}
 * 0 & 0 & 0 & -L & 0 \\
 * 0 & 0 & 0 & 0 & -C \\
 * 0 & 0 & 0 & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccccc}
 * 1 & -R_{ste}-R & -1 & 0 & 0 \\
 * 0 & 1 & 0 & -1 & 0 \\
 * 1 & -R_{ste}-R & 0 & 0 & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{D}^{e}=\left[\begin{array}{ccccc}
 * 0 & -R_{ste} & 0 & 0 & 0 \\
 * 0 & 0 & 0 & 0 & 0 \\
 * 0 & -R_{ste} & 0 & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * with the stenosis resistance \f$ R_{ste}=K_{t} \frac{\rho}{2
 * A_{o}^{2}}\left(\frac{A_{o}}{A_{s}}-1\right)^{2}|Q_{in}^{e}| \f$. The
 * constant part of the equation is summarized in \ref
 * Parameters::stenosis_coefficient. \f$R\f$, \f$C\f$, and \f$L\f$ refer to
 * Poisieuille resistance, capacitance and inductance, respectively.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Poiseuille resistance
 * * `1` Capacitance
 * * `2` Inductance
 * * `3` Stenosis coefficient
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class BloodVessel : public Block<T> {
 public:
  // Inherit constructors
  using Block<T>::Block;

  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    RESISTANCE = 0,
    CAPACITANCE = 1,
    INDUCTANCE = 2,
    STENOSIS_COEFFICIENT = 3,
  };

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
  void update_solution(ALGEBRA::SparseSystem<T> &system,
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
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &dy);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  std::map<std::string, int> num_triplets = {
      {"F", 5},
      {"E", 3},
      {"D", 2},
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
void BloodVessel<T>::setup_dofs(DOFHandler &dofhandler) {
  Block<T>::setup_dofs_(dofhandler, 2, {});
}

template <typename T>
void BloodVessel<T>::update_constant(ALGEBRA::SparseSystem<T> &system,
                                     std::vector<T> &parameters) {
  // Get parameters
  T capacitance = parameters[this->global_param_ids[ParamId::CAPACITANCE]];
  T inductance = parameters[this->global_param_ids[ParamId::INDUCTANCE]];

  // Set element contributions
  system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[3]) =
      -inductance;
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) =
      -capacitance;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[3]) = -1.0;
}

template <typename T>
void BloodVessel<T>::update_solution(ALGEBRA::SparseSystem<T> &system,
                                     std::vector<T> &parameters,
                                     Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
                                     Eigen::Matrix<T, Eigen::Dynamic, 1> &dy) {
  // Get parameters
  T resistance = parameters[this->global_param_ids[ParamId::RESISTANCE]];
  T capacitance = parameters[this->global_param_ids[ParamId::CAPACITANCE]];
  T stenosis_coeff =
      parameters[this->global_param_ids[ParamId::STENOSIS_COEFFICIENT]];
  T q_in = y[this->global_var_ids[1]];
  T dq_in = dy[this->global_var_ids[1]];
  T stenosis_resistance = stenosis_coeff * fabs(q_in);

  // Set element contributions
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
      capacitance * (resistance + 2.0 * stenosis_resistance);
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
      -resistance - stenosis_resistance;
  system.D.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
      -stenosis_resistance;

  T sgn_q_in = (0.0 < q_in) - (q_in < 0.0);
  system.D.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
      2.0 * capacitance * stenosis_coeff * sgn_q_in * dq_in;
}

template <typename T>
void BloodVessel<T>::update_gradient(
    Eigen::SparseMatrix<T> &jacobian,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &residual,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &alpha,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &dy) {
  T y0 = y[this->global_var_ids[0]];
  T y1 = y[this->global_var_ids[1]];
  T y2 = y[this->global_var_ids[2]];
  T y3 = y[this->global_var_ids[3]];

  T dy0 = dy[this->global_var_ids[1]];
  T dy1 = dy[this->global_var_ids[1]];
  T dy3 = dy[this->global_var_ids[3]];

  T resistance = alpha[this->global_param_ids[ParamId::RESISTANCE]];
  T capacitance = alpha[this->global_param_ids[ParamId::CAPACITANCE]];
  T inductance = alpha[this->global_param_ids[ParamId::INDUCTANCE]];
  T stenosis_coeff = 0.0;
  if (this->global_eqn_ids.size() > 3) {
    stenosis_coeff =
        alpha[this->global_param_ids[ParamId::STENOSIS_COEFFICIENT]];
  }
  T stenosis_resistance = stenosis_coeff * fabs(y1);

  jacobian.coeffRef(this->global_eqn_ids[0], this->global_param_ids[0]) = -y1;
  jacobian.coeffRef(this->global_eqn_ids[0], this->global_param_ids[2]) = -dy3;

  if (this->global_eqn_ids.size() > 3) {
    jacobian.coeffRef(this->global_eqn_ids[0], this->global_param_ids[3]) =
        -fabs(y1) * y1;
  }

  jacobian.coeffRef(this->global_eqn_ids[1], this->global_param_ids[0]) =
      capacitance * y1 * dy1;
  jacobian.coeffRef(this->global_eqn_ids[1], this->global_param_ids[1]) =
      -y0 + (resistance + 2 * stenosis_resistance) * dy1;

  if (this->global_eqn_ids.size() > 3) {
    jacobian.coeffRef(this->global_eqn_ids[1], this->global_param_ids[3]) =
        2.0 * capacitance * fabs(y1) * dy1;
  }

  residual(this->global_eqn_ids[0]) =
      y0 - (resistance + stenosis_resistance) * y1 - y2 - inductance * dy3;
  residual(this->global_eqn_ids[1]) =
      y1 - y3 - capacitance * dy0 +
      capacitance * (resistance + 2 * stenosis_resistance) * y1 * dy1;
}

template <typename T>
std::map<std::string, int> BloodVessel<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_BLOODVESSEL_HPP_
