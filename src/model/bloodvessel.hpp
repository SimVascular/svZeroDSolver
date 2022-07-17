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

#include "../algebra/densesystem.hpp"
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
 * to [R, l=$R_{ste}$, -] (5,0)
 * to [L, l=$L$, *-*] (7,0)
 * node[anchor=south]{$P_{out}$}
 * (5,0) to [C, l=$C$, *-] (5,-1.5)
 * node[ground]{$P_C$};
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
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class BloodVessel : public Block<T> {
 public:
  /**
   * @brief Parameters of the element.
   *
   * Struct containing all constant and/or time-dependent parameters of the
   * element.
   */
  struct Parameters : public Block<T>::Parameters {
    T R;                     ///< Poiseuille resistance
    T C;                     ///< Capacitance
    T L;                     ///< Inductance
    T stenosis_coefficient;  ///< Stenosis Coefficient
  };

  /**
   * @brief Construct a new Blood Vessel object
   *
   * @param R Poiseuille resistance
   * @param C Capacitance
   * @param L Inductance
   * @param stenosis_coefficient Stenosis Coefficient
   * @param name Name
   */
  BloodVessel(T R, T C, T L, T stenosis_coefficient, std::string name);

  /**
   * @brief Destroy the Blood Vessel object
   *
   */
  ~BloodVessel();

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
   * @brief Update the constant contributions of the element in a dense system
   *
   * @param system System to update contributions at
   */
  void update_constant(ALGEBRA::DenseSystem<T> &system);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * dense system
   *
   * @param system System to update contributions at
   * @param y Current solution
   */
  void update_solution(ALGEBRA::DenseSystem<T> &system,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

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
      {"F", 10},
      {"E", 2},
      {"D", 2},
  };

 private:
  Parameters params;
};

template <typename T>
BloodVessel<T>::BloodVessel(T R, T C, T L, T stenosis_coefficient,
                            std::string name)
    : Block<T>(name) {
  this->name = name;
  this->params.R = R;
  this->params.C = C;
  this->params.L = L;
  this->params.stenosis_coefficient = stenosis_coefficient;
}

template <typename T>
BloodVessel<T>::~BloodVessel() {}

template <typename T>
void BloodVessel<T>::setup_dofs(DOFHandler &dofhandler) {
  Block<T>::setup_dofs_(dofhandler, 3, 1);
}

template <typename T>
void BloodVessel<T>::update_constant(ALGEBRA::DenseSystem<T> &system) {
  system.E(this->global_eqn_ids[0], this->global_var_ids[3]) = -params.L;
  system.E(this->global_eqn_ids[1], this->global_var_ids[4]) = -params.C;

  system.F(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
  system.F(this->global_eqn_ids[0], this->global_var_ids[1]) = -params.R;
  system.F(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;

  system.F(this->global_eqn_ids[1], this->global_var_ids[1]) = 1.0;
  system.F(this->global_eqn_ids[1], this->global_var_ids[3]) = -1.0;

  system.F(this->global_eqn_ids[2], this->global_var_ids[0]) = 1.0;
  system.F(this->global_eqn_ids[2], this->global_var_ids[1]) = -params.R;
  system.F(this->global_eqn_ids[2], this->global_var_ids[4]) = -1.0;
}

template <typename T>
void BloodVessel<T>::update_solution(ALGEBRA::DenseSystem<T> &system,
                                     Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  T q_in = fabs(y[this->inlet_nodes[0]->flow_dof]);
  T fac1 = -params.stenosis_coefficient * q_in;
  T fac2 = fac1 - params.R;
  system.F(this->global_eqn_ids[0], this->global_var_ids[1]) = fac2;
  system.F(this->global_eqn_ids[2], this->global_var_ids[1]) = fac2;
  system.D(this->global_eqn_ids[0], this->global_var_ids[1]) = fac1;
  system.D(this->global_eqn_ids[2], this->global_var_ids[1]) = fac1;
}

template <typename T>
void BloodVessel<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[3]) =
      -params.L;
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[4]) =
      -params.C;

  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
      -params.R;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;

  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[3]) = -1.0;

  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[0]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[1]) =
      -params.R;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[4]) = -1.0;
}

template <typename T>
void BloodVessel<T>::update_solution(ALGEBRA::SparseSystem<T> &system,
                                     Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  T q_in = fabs(y[this->inlet_nodes[0]->flow_dof]);
  T fac1 = -params.stenosis_coefficient * q_in;
  T fac2 = fac1 - params.R;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = fac2;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[1]) = fac2;
  system.D.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = fac1;
  system.D.coeffRef(this->global_eqn_ids[2], this->global_var_ids[1]) = fac1;
}
}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_BLOODVESSEL_HPP_