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
 * @file openloopcoronarybc.hpp
 * @brief MODEL::OpenLoopCoronaryBC source file
 */
#ifndef SVZERODSOLVER_MODEL_OPENLOOPCORONARYBC_HPP_
#define SVZERODSOLVER_MODEL_OPENLOOPCORONARYBC_HPP_

#include "../algebra/sparsesystem.hpp"
#include "block.hpp"
#include "parameter.hpp"

namespace MODEL {

/**
 * @brief Open loop coronary boundary condition based on \cite kim_coronary.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_a$, *-] (3,0)
 * to [R, l=$R_{am}$, -] (5,0)
 * to [R, l=$R_v$, *-*] (7,0)
 * node[anchor=south]{$P_{v}$}
 * (5,0) to [C, l=$C_{im} \;V_{im}$, -*] (5,-1.5)
 * node[left]{$P_{im}$}
 * (3,0) to [C, l=$C_a$, -*] (3,-1.5)
 * node[left]{$P_a$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * C_{i m} R_{v} Q^{e}-V_{i m}^{e}-C_{i m} P_{i m}+C_{i m} P_{v}-C_{i m} R_{v}
 * \frac{d V_{i m}^{e}}{d t}-C_{a} C_{i m} R_{v} \frac{d P^{e}}{d t}+R_{a} C_{a}
 * C_{i m} R_{v} \frac{d Q^{e}}{d t}+C_{a} C_{i m} R_{v} \frac{d P_{a}^{e}}{d
 * t}=0 \f]
 *
 * \f[
 * C_{i m} R_v P^{e}-C_{i m} R_{v} R_{a} Q^{e}-R_{v} V_{i m}^{e}-C_{i m} R_{v}
 * P_{i m}-C_{i m} R_{v} R_{a m} \frac{d V_{i m}^{e}}{d t}-R_{a m} V_{i
 * m}^{e}-C_{i m} R_{a m} P_{i m}+R_{a m} C_{i m} P_{v}=0 \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lll}P^{e} & Q^{e} & V_{i
 * m}^{e}\end{array}\right]^{T}, \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{ccc}-C_{a} C_{i m} R_{v} & R_{a} C_{a}
 * C_{i m} R_{v} & -C_{i m} R_{v} \\ 0 & 0 & -C_{i m} R_{v} R_{a
 * m}\end{array}\right] \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccc}0 & C_{i m} R_{v} & -1 \\C_{i m} R_{v}
 * & -C_{i m} R_{v} R_{a} & -\left(R_{v}+R_{a m}\right)\end{array}\right] \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}C_{i m}\left(-P_{i m}+P_{v}\right)+C_{a}
 * C_{i m} R_{v} \frac{d P_{a}}{d t} \\-C_{i m}\left(R_{v}+R_{a m}\right) P_{i
 * m}+R_{a m} C_{i m} P_{v}\end{array}\right] \f]
 *
 * Assume \f$P_a=0\f$.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Ra
 * * `1` Ram
 * * `2` Rv
 * * `3` Ca
 * * `4` Cim
 * * `5` Pim
 * * `6` Pv
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class OpenLoopCoronaryBC : public Block<T> {
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
      {"F", 5},
      {"E", 4},
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
void OpenLoopCoronaryBC<T>::setup_dofs(DOFHandler &dofhandler) {
  Block<T>::setup_dofs_(dofhandler, 2, {"volume_im"});
}

template <typename T>
void OpenLoopCoronaryBC<T>::update_constant(ALGEBRA::SparseSystem<T> &system,
                                            std::vector<T> &parameters) {
  T Ra = parameters[this->global_param_ids[0]];
  T Ram = parameters[this->global_param_ids[1]];
  T Rv = parameters[this->global_param_ids[2]];
  T Ca = parameters[this->global_param_ids[3]];
  T Cim = parameters[this->global_param_ids[4]];
  if (this->steady) {
    // Different assmembly for steady block to avoid singular system
    // and solve for the internal variable V_im inherently
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = -Cim;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
        Cim * (Ra + Ram);
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = 1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) = -1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
        Ra + Ram + Rv;
  } else {
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
        Cim * Rv;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) =
        Cim * Rv;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
        -Cim * Rv * Ra;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) =
        -(Rv + Ram);

    system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) =
        -Ca * Cim * Rv;
    system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
        Ra * Ca * Cim * Rv;
    system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) =
        -Cim * Rv;
    system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) =
        -Cim * Rv * Ram;
  }
}

template <typename T>
void OpenLoopCoronaryBC<T>::update_time(ALGEBRA::SparseSystem<T> &system,
                                        std::vector<T> &parameters) {
  T Ram = parameters[this->global_param_ids[1]];
  T Rv = parameters[this->global_param_ids[2]];
  T Cim = parameters[this->global_param_ids[4]];
  T Pim = parameters[this->global_param_ids[5]];
  T Pv = parameters[this->global_param_ids[6]];
  if (this->steady) {
    system.C(this->global_eqn_ids[0]) = -Cim * Pim;
    system.C(this->global_eqn_ids[1]) = Pv;
  } else {
    system.C(this->global_eqn_ids[0]) = Cim * (-Pim + Pv);
    system.C(this->global_eqn_ids[1]) =
        -Cim * (Rv + Ram) * Pim + Ram * Cim * Pv;
  }
}

template <typename T>
std::map<std::string, int> OpenLoopCoronaryBC<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_OPENLOOPCORONARYBC_HPP_
