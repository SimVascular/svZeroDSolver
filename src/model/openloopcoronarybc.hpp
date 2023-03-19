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
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class OpenLoopCoronaryBC : public Block<T> {
 public:
  /**
   * @brief Parameters of the element.
   *
   * Struct containing all constant and/or time-dependent parameters of the
   * element.
   */
  struct Parameters : public Block<T>::Parameters {
    T Ra;
    T Ram;
    T Rv;
    T Ca;
    T Cim;
    Parameter<T> Pim;  ///<
    Parameter<T> Pv;
  };

  /**
   * @brief Construct a new OpenLoopCoronaryBC object
   *
   * @param P Time dependent pressure
   * @param name Name
   */
  OpenLoopCoronaryBC(T Ra, T Ram, T Rv, T Ca, T Cim, Parameter<T> Pim,
                     Parameter<T> Pv, std::string name);

  /**
   * @brief Destroy the OpenLoopCoronaryBC object
   *
   */
  ~OpenLoopCoronaryBC();

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

  /**
   * @brief Convert the block to a steady behavior
   *
   * Converts the prescribed pressure to the constant mean of itself
   *
   */
  void to_steady();

  /**
   * @brief Convert the block to an unsteady behavior
   *
   */
  void to_unsteady();

 private:
  Parameters params;
  bool issteady = false;
};

template <typename T>
OpenLoopCoronaryBC<T>::OpenLoopCoronaryBC(T Ra, T Ram, T Rv, T Ca, T Cim,
                                          Parameter<T> Pim, Parameter<T> Pv,
                                          std::string name)
    : Block<T>(name) {
  this->name = name;
  this->params.Ra = Ra;
  this->params.Ram = Ram;
  this->params.Rv = Rv;
  this->params.Ca = Ca;
  this->params.Cim = Cim;
  this->params.Pim = Pim;
  this->params.Pv = Pv;
}

template <typename T>
OpenLoopCoronaryBC<T>::~OpenLoopCoronaryBC() {}

template <typename T>
void OpenLoopCoronaryBC<T>::setup_dofs(DOFHandler &dofhandler) {
  Block<T>::setup_dofs_(dofhandler, 2, {"volume_im"});
}

template <typename T>
void OpenLoopCoronaryBC<T>::update_block_params(std::vector<T> new_params) {
  this->params.Ra = new_params[0];
  this->params.Ram = new_params[1];
  this->params.Rv = new_params[2];
  this->params.Ca = new_params[3];
  this->params.Cim = new_params[4];
  // this->params.Pim = params[5];
  // this->params.Pv = params[6];
  int num_time_pts = (int)new_params[5];
  std::vector<T> t_new;
  std::vector<T> Pv_new;
  std::vector<T> Pim_new;
  for (int i = 0; i < num_time_pts; i++) {
    t_new.push_back(new_params[1 + i]);
    Pv_new.push_back(new_params[1 + num_time_pts + i]);
    Pim_new.push_back(new_params[1 + 2 * num_time_pts + i]);
  }
  this->params.Pv.update_params(t_new, Pv_new);
  this->params.Pim.update_params(t_new, Pim_new);
}

template <typename T>
void OpenLoopCoronaryBC<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  if (issteady) {
    // Different assmembly for steady block to avoid singular system
    // and solve for the internal variable V_im inherently
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) =
        -params.Cim;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
        params.Cim * (params.Ra + params.Ram);
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = 1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) = -1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
        params.Ra + params.Ram + params.Rv;
  } else {
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
        params.Cim * params.Rv;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) =
        params.Cim * params.Rv;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
        -params.Cim * params.Rv * params.Ra;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) =
        -(params.Rv + params.Ram);

    system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) =
        -params.Ca * params.Cim * params.Rv;
    system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
        params.Ra * params.Ca * params.Cim * params.Rv;
    system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) =
        -params.Cim * params.Rv;
    system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) =
        -params.Cim * params.Rv * params.Ram;
  }
}

template <typename T>
void OpenLoopCoronaryBC<T>::update_time(ALGEBRA::SparseSystem<T> &system,
                                        T time) {
  T Pim = params.Pim.get(time);
  T Pv = params.Pv.get(time);

  if (issteady) {
    system.C(this->global_eqn_ids[0]) = -params.Cim * Pim;
    system.C(this->global_eqn_ids[1]) = Pv;
  } else {
    system.C(this->global_eqn_ids[0]) = params.Cim * (-Pim + Pv);
    system.C(this->global_eqn_ids[1]) =
        -params.Cim * (params.Rv + params.Ram) * Pim +
        params.Ram * params.Cim * Pv;
  }
}

template <typename T>
void OpenLoopCoronaryBC<T>::to_steady() {
  params.Pim.to_steady();
  params.Pv.to_steady();
  issteady = true;
}

template <typename T>
void OpenLoopCoronaryBC<T>::to_unsteady() {
  params.Pim.to_unsteady();
  params.Pv.to_unsteady();
  issteady = false;
}

template <typename T>
std::map<std::string, int> OpenLoopCoronaryBC<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_OPENLOOPCORONARYBC_HPP_
