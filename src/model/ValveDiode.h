// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ValveDiode.h
 * @brief model::ValveDiode source file
 */
#ifndef SVZERODSOLVER_MODEL_VALVEDIODE_HPP_
#define SVZERODSOLVER_MODEL_VALVEDIODE_HPP_

#include <math.h>

#include "Block.h"
#include "SparseSystem.h"
#include "debug.h"

/**
 * @brief Valve (regularized diode) block.
 *
 * Models a check valve as a resistor with two slopes, \f$R_{min}\f$
 * forward and \f$R_{max}\f$ in reverse,
 *
 * \f[
 * Q = \frac{\max(\Delta p, 0)}{R_{min}} + \frac{\min(\Delta p, 0)}
 * {R_{max}}, \qquad \Delta p = P_{in} - P_{out}.
 * \f]
 *
 * A stenosis raises \f$R_{min}\f$ and a regurgitant valve lowers
 * \f$R_{max}\f$. Leaving \f$R_{max}\f$ out gives the ideal check valve
 * \f$R_{min} Q = \max(\Delta p, 0)\f$, which is the complementarity
 * condition \f$Q \geq 0\f$, \f$R_{min} Q - \Delta p \geq 0\f$,
 * \f$Q\,(R_{min} Q - \Delta p) = 0\f$.
 *
 * Replacing the maximum by the smoothed
 * \f$\max_{\varepsilon}(x) = \frac{1}{2}(x + \sqrt{x^2 +
 * 4\varepsilon^2})\f$ and writing \f$\beta = R_{min}/R_{max} \in
 * [0,1)\f$ for the reverse-to-forward conductance ratio, the governing
 * equation becomes
 *
 * \f[
 * R_{min}Q_{in} - \frac{1}{2}(1+\beta)\,\Delta p
 * - \frac{1}{2}(1-\beta)\sqrt{\Delta p^2 + 4\varepsilon^2} = 0.
 * \f]
 *
 * \f$\varepsilon\f$ is a pressure and sets the width of the transition
 * between open and closed. At \f$\varepsilon = 0\f$ the exact two-slope
 * valve is recovered and Newton becomes a semi-smooth (active-set)
 * iteration; for \f$\varepsilon > 0\f$ the equation is smooth and the
 * valve leaks at most \f$\varepsilon / R_{min}\f$ at \f$\Delta p =
 * 0\f$, decaying like \f$\varepsilon^2 / (R_{min}|\Delta p|)\f$ once
 * firmly closed.
 *
 * Compared with \ref ValveTanh, the derivative with respect to
 * \f$\Delta p\f$ stays between \f$-1\f$ and \f$-\beta\f$ whatever
 * \f$R_{max}/R_{min}\f$ is, instead of growing with
 * \f$k\,Q\,(R_{max}-R_{min})\f$ at the transition; there is no steepness
 * whose units have to be tracked across unit systems; and the solution
 * converges to the ideal valve as \f$\varepsilon \to 0\f$ instead of
 * changing qualitatively. Since \f$\beta\f$ enters linearly, the equation
 * stays piecewise linear for every value of it, so regurgitation makes the
 * problem less non-linear rather than more. At \f$\beta = 1\f$ the valve
 * would be a plain resistor, which the block rejects.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [D, l=$R_{min}$, *-*] (3,0)
 * node[anchor=south]{$P_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * R_{min}Q_{in} - \frac{1}{2}(1+\beta)\left(P_{in}-P_{out}\right)
 * - \frac{1}{2}(1-\beta)\sqrt{(P_{in}-P_{out})^2 + 4\varepsilon^2} = 0
 * \f]
 *
 * \f[
 * Q_{in}-Q_{out}=0
 * \f]
 *
 * \f[
 * \text{valve\_status} = \frac{1}{2}\left(1 -
 * \frac{\Delta p}{\sqrt{\Delta p^2 + 4\varepsilon^2}}\right)
 * \f]
 *
 * where valve_status is 0 for a fully open and 1 for a fully closed valve,
 * matching the convention of \ref ValveTanh.
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{llll}P_{in} & Q_{in} &
 * P_{out} & Q_{out} \end{array}\right]^{T} \f]
 *
 * \f[
 * \mathbf{E}^{e}=\mathbf{0}
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{cccc}
 * -(1+\beta)/2 & R_{min} &  (1+\beta)/2 & 0\\
 * 0            & 1       &  0            & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * -\frac{1}{2}(1-\beta)\sqrt{\Delta p^2 + 4\varepsilon^2} \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \left(\frac{\partial\mathbf{c}}{\partial\mathbf{y}}\right)^{e} =
 * \left[\begin{array}{cccc}
 * -d & 0 &  d & 0\\
 * 0  & 0 &  0 & 0
 * \end{array}\right] \f]
 * with
 * \f[
 * s = \frac{1}{2}\left(1 +
 * \frac{\Delta p}{\sqrt{\Delta p^2 + 4\varepsilon^2}}\right)
 * \in [0, 1], \qquad d = (1-\beta)\left(s - \frac{1}{2}\right).
 * \f]
 *
 * The valve status is a post-processing variable, see
 * DOFHandler::register_output_variable: it has no equation and is computed
 * from the converged state in \ref update_outputs.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Rmin: Resistance of the open valve, raised by a stenosis (must be
 *   positive)
 * * `1` Rmax: Resistance of the closed valve, lowered by a regurgitant
 *   valve. Optional and must exceed `Rmin`; left out, the valve does not
 *   leak.
 * * `2` eps: Smoothing pressure, optional, defaults to 0 (exact valve)
 * * `3` upstream_block: Name of block connected upstream
 * * `4` downstream_block: Name of block connected downstream
 *
 * ### Usage in json configuration file
 *
 *     "valves": [
 *         {
 *             "type": "ValveDiode",
 *             "name": "valve",
 *             "params": {
 *                 "Rmin": 0.005,
 *                 "upstream_block": "upstream_vessel",
 *                 "downstream_block": "downstream_vessel"
 *             }
 *         }
 *     ]
 *
 * A stenotic valve raises `Rmin`; a regurgitant one adds a finite `Rmax`.
 *
 * ### Internal variables
 *
 * This block has no internal variables. It reports one post-processing
 * variable, `valve_status`.
 *
 */
class ValveDiode : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    RMIN = 0,
    RMAX = 1,
    EPS = 2,
  };

  /**
   * @brief Construct a new ValveDiode object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ValveDiode(int id, Model* model)
      : Block(id, model, BlockType::valve_diode, BlockClass::valve,
              {{"Rmin", InputParameter()},
               {"Rmax", InputParameter(true, false, true, 0.0)},
               {"eps", InputParameter(true, false, true, 0.0)},
               {"upstream_block", InputParameter(false, false, false)},
               {"downstream_block", InputParameter(false, false, false)}}) {}

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler& dofhandler) override;

  /**
   * @brief Check the resistances once, after the parameters are loaded
   */
  void setup_model_dependent_params() override;

  /**
   * @brief Update the constant contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem& system,
                       std::vector<double>& parameters) override;

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(
      SparseSystem& system, std::vector<double>& parameters,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) override;

  /**
   * @brief Compute the valve status from the converged state
   *
   * @param y Accepted state vector
   * @param ydot Accepted derivative of the state vector
   */
  void update_outputs(
      Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& ydot) override;

 private:
  /**
   * @brief Reverse-to-forward conductance ratio Rmin / Rmax
   *
   * Zero for a valve that does not leak, i.e. a non-positive Rmax.
   *
   * @param parameters Parameters of the model
   * @return Conductance ratio in [0, 1)
   */
  double get_beta(const std::vector<double>& parameters) const;

 public:
  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{5, 0, 2};
};

#endif  // SVZERODSOLVER_MODEL_VALVEDIODE_HPP_
