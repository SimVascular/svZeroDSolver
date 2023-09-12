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

#include "../algebra/SparseSystem.h"
#include "Block.h"

namespace zd_model {

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
 * (5,0) to [L, l=$L$, -*] (7,0)
 * node[anchor=south]{$P_{out}$}
 * (5,0) to [C, l=$C$, -] (5,-1.5)
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
 * Q_{i n}^{e}-Q_{o u t}^{e}-C \frac{d P_{in}^{e}}{d t}+C(R+2R_{ste})\frac{d
 * Q_{in}^{e}}{d t}=0 \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{llll}P_{i n}^{e} & Q_{in}^{e} &
 * P_{out}^{e} & Q_{out}^{e}\end{array}\right]^{T} \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{cccc}
 * 1 & -R_{ste}-R & -1 & 0 \\
 * 0 & 1 & 0 & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{cccc}
 * 0 & 0 & 0 & -L \\
 * -C & C(R+2R_{ste}) & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{D}^{e}=\left[\begin{array}{cccc}
 * 0 & -R_{ste} & 0 & 0 \\
 * 0 & 2CK_{ste} sgn(Q_{in}^{e}) \dot{Q}_{in}^{e} & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * with the stenosis resistance \f$ R_{ste}=K_{t} \frac{\rho}{2
 * A_{o}^{2}}\left(\frac{A_{o}}{A_{s}}-1\right)^{2}|Q_{in}^{e}| \f$. The
 * constant part of the equation is summarized in \ref
 * Parameters::stenosis_coefficient. \f$R\f$, \f$C\f$, and \f$L\f$ refer to
 * Poisieuille resistance, capacitance and inductance, respectively.
 *
 * ### Gradient
 *
 * Gradient of the equations with respect to the parameters:
 *
 * \f[
 * \mathbf{J}^{e} = \left[\begin{array}{cccc}
 * -y_2 & 0 & -\dot{y}_4 & -|y_2|y_2 \\
 * C\dot{y}_2 & (-\dot{y}_1+(R+2R_{ste})\dot{y}_2) & 0 & 2C|y_2|\dot{y}_2
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{r}^{e} = \left[\begin{array}{c}
 * y_1-(R+R_{ste})y_2-y_3-L\dot{y}_4 \\
 * y_2 - y_4 - C\dot{y}_1 + C(R+2R_{ste}) \dot{y}_2
 * \end{array}\right]
 * \f]
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
class BloodVessel : public Block {
 public:
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

  explicit BloodVessel(int id, const std::vector<int> &param_ids, Model *model) : Block(id, param_ids, model){};


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
  void update_constant(algebra::SparseSystem& system, std::vector<double> &parameters);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(algebra::SparseSystem& system, std::vector<double> &parameters,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

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
  void update_gradient(Eigen::SparseMatrix<double>& jacobian,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &residual,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &alpha,
                       std::vector<double> &y, std::vector<double> &dy);

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

}  

#endif  // SVZERODSOLVER_MODEL_BLOODVESSEL_HPP_