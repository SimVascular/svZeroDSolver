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
 * @file closedloopRCRbc.hpp
 * @brief MODEL::ClosedLoopRCRBC source file
 */
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPRCRBC_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPRCRBC_HPP_

#include "../algebra/sparsesystem.hpp"
#include "block.hpp"

namespace MODEL {
/**
 * @brief Closed-loop RCR boundary condition.
 *
 * Models the mechanical behavior of a Windkessel boundary condition that is
 * connected to other blocks on both sides.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_p$, *-] (3,0)
 * to [R, l=$R_d$, *-*] (5,0)
 * node[anchor=south]{$P_{out}$}
 * (3,0) to [C, l=$C$, *-] (3,-1.5)
 * node[ground]{$P_{C}$};
 * \draw [-latex] (5.2,0) -- (6.0,0) node[right] {$Q_{out}$} ;
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * C \frac{d P_c}{dt} + Q_{out} -Q_{in} = 0
 * \f]
 *
 * \f[
 * P_{in}-P_{c}-R_{p} Q_{in}=0
 * \f]
 *
 * \f[
 * P_{c} - P_{out} - R_{d} Q_{out}=0
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}=\left[\begin{array}{lllll}P_{in} & Q_{in} & P_{out} & Q_{out} &
 * P_{c}\end{array}\right]^{T} \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{ccccc}
 * 0 & 0 & 0 & 0 & C \\
 * 0 & 0 & 0 & 0 & 0 \\
 * 0 & 0 & 0 & 0 & 0 \\
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccccc}
 * 0 & -1 & 1 & 0 & 0 \\
 * 1 & -R_p & 0 & 0 & -1 \\
 * 0 & 0 & -1 & -R_d & +1 \\
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * 0 \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class ClosedLoopRCRBC : public Block<T> {
 public:
  /**
   * @brief Parameters of the element.
   *
   * Struct containing all constant and/or time-dependent parameters of the
   * element.
   */
  struct Parameters : public Block<T>::Parameters {
    T Rp;  ///< Proximal resistance
    T C;   ///< Capacitance
    T Rd;  ///< Distal restistance
  };

  /**
   * @brief Construct a new ClosedLoopRCRBC object
   *
   * @param Rp Proximal resistance
   * @param C Capacitance
   * @param Rd Distal resistance
   * @param Pd Distal pressure
   * @param name Name
   */
  ClosedLoopRCRBC(T Rp, T C, T Rd, bool closed_loop_outlet, std::string name);

  /**
   * @brief Destroy the ClosedLoopRCRBC object
   *
   */
  ~ClosedLoopRCRBC();

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
      {"F", 8},
      {"E", 1},
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
   * Set the capacitance to 0.
   */
  void to_steady();

  /**
   * @brief Convert the block to a steady behavior
   *
   * Set the capacitance to original value.
   */
  void to_unsteady();

 private:
  Parameters params;
  T c_cache;
  bool closed_loop_outlet =
      false;  ///< Is this block connected to a closed-loop model?
};

template <typename T>
ClosedLoopRCRBC<T>::ClosedLoopRCRBC(T Rp, T C, T Rd, bool closed_loop_outlet,
                                    std::string name)
    : Block<T>(name) {
  this->name = name;
  this->params.Rp = Rp;
  this->params.C = C;
  this->params.Rd = Rd;
  this->closed_loop_outlet = closed_loop_outlet;
}

template <typename T>
ClosedLoopRCRBC<T>::~ClosedLoopRCRBC() {}

template <typename T>
void ClosedLoopRCRBC<T>::setup_dofs(DOFHandler &dofhandler) {
  // Block<T>::setup_dofs_(dofhandler, 3, 1);
  Block<T>::setup_dofs_(dofhandler, 3, {"P_c"});
}

template <typename T>
void ClosedLoopRCRBC<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[3]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[4]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[2]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[4]) = 1.0;

  // Below values can be unsteady if needed (not currently implemented)
  system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[4]) =
      params.C;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
      -params.Rp;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[3]) =
      -params.Rd;
}

template <typename T>
void ClosedLoopRCRBC<T>::to_steady() {
  c_cache = params.C;
  params.C = 0.0;
}

template <typename T>
void ClosedLoopRCRBC<T>::to_unsteady() {
  params.C = c_cache;
}

template <typename T>
std::map<std::string, int> ClosedLoopRCRBC<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_CLOSEDLOOPRCRBCBC_HPP_
