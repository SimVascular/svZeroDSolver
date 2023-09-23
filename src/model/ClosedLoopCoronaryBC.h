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
 * @file closedloopcoronarybc.hpp
 * @brief MODEL::ClosedLoopCoronaryBC source file
 */
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYBC_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYBC_HPP_

#include "Block.h"
#include "ClosedLoopHeartPulmonary.h"
#include "SparseSystem.h"

namespace zd_model {

enum class Side { LEFT, RIGHT, NONE };

/**
 * @brief Closed loop coronary boundary condition (connected to other blocks on
 * both sides).
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_a$, *-] (3,0)
 * to [R, l=$R_{am}$, -] (5,0)
 * to [R, l=$R_v$, *-*] (7,0)
 * node[anchor=south]{$P_{out}$}
 * (5,0) to [C, l=$C_{im} \;V_{im}$, -*] (5,-1.5)
 * node[left]{$P_{im}$}
 * (3,0) to [C, l=$C_a$, -*] (3,-1.5)
 * node[left]{$P_a$};
 * \draw [-latex] (7.2,0) -- (8.0,0) node[right] {$Q_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_{out} - P_{in} + (R_{am}+R_a)Q_{in} + R_v Q_{out} + R_{am} C_a
 * \frac{dP_a}{dt} - R_{am} C_a \frac{dP_{in}}{dt} + R_{am} R_a C_a
 * \frac{dQ_{in}}{dt} = 0 \f]
 *
 * \f[
 * Q_{in} - Q_{out} + C_a \frac{dP_a}{dt} - C_a \frac{dP_{in}}{dt} + C_a R_a
 * \frac{dQ_{in}}{dt} - \frac{dV_{im}}{dt} = 0 \f]
 *
 * \f[
 * C_{im} P_{out} + C_{im} R_v Q_{out} - C_{im} P_{im} - V_{im} = 0
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lllll}P^{in} & Q^{in} & P_{out} & Q_{out}
 * & V_{im}^{e}\end{array}\right]^{T}, \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{ccccc}
 * -R_{am} C_{a} & R_{am} R_{a} C_{a} & 0 & 0 & 0 \\
 * -C_{a} & R_{a} C_{a} & 0 & 0 & -1 \\
 * 0 & 0 & 0 & 0 & 0 \\
 * \end{array}\right] \f]
 *
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccccc}
 * -1 & R_{am} + R_{a} & 1 & R_v & 0 \\
 * 0 & 1 & 0 & -1 & 0 \\
 * 0 & 0 & C_{im} & C_{im} R_v & -1 \\
 * \end{array}\right] \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * C_{a} R_{am} \frac{d P_{a}}{d t} \\
 * C_{a}\frac{d P_{a}}{d t} \\
 * -C_{im} P_{im}
 * \end{array}\right] \f]
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
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 * @tparam side Side of the block (e.g. `Side::LEFT`, `Side::RIGHT`)
 */
class ClosedLoopCoronaryBC : public Block {
 public:
  explicit ClosedLoopCoronaryBC(int id, const std::vector<int> &param_ids,
                                Model *model, Side side)
      : Block(id, param_ids, model), side{side} {};

  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId { RA = 0, RAM = 1, RV = 2, CA = 3, CIM = 4 };

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
   * @brief Setup parameters that depend on the model
   *
   */
  void setup_model_dependent_params();

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(algebra::SparseSystem &system,
                       std::vector<double> &parameters);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(algebra::SparseSystem &system,
                       std::vector<double> &parameters,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  std::map<std::string, int> num_triplets = {
      {"F", 9},
      {"E", 5},
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
  int ventricle_var_id;  // Variable index of either left or right ventricle
  int im_param_id;       // Index of parameter Im
  Side side{Side::NONE};
};

}  // namespace zd_model

#endif  // SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYBC_HPP_
