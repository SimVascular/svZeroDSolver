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
 * Models an ideal check valve with open resistance \f$R\f$ as the
 * complementarity (diode) condition
 *
 * \f[
 * Q \geq 0, \quad RQ - \Delta p \geq 0, \quad Q\,(RQ - \Delta p) = 0,
 * \qquad \Delta p = P_{in} - P_{out},
 * \f]
 *
 * which is equivalent to the single equation \f$ RQ = \max(\Delta p, 0)\f$.
 * The valve is open when \f$\Delta p > 0\f$ (then \f$\Delta p = RQ\f$) and
 * closed otherwise (then \f$Q = 0\f$).
 *
 * The non-smooth maximum is replaced by the standard smoothed maximum
 *
 * \f[
 * \max\nolimits_{\varepsilon}(x) =
 * \frac{1}{2}\left(x + \sqrt{x^2 + 4\varepsilon^2}\right),
 * \f]
 *
 * so that the governing equation reads
 *
 * \f[
 * RQ_{in} - \frac{1}{2}\left(\Delta p +
 * \sqrt{\Delta p^2 + 4\varepsilon^2}\right) = 0.
 * \f]
 *
 * \f$\varepsilon\f$ is a pressure and sets the width of the transition
 * between the open and the closed state. For \f$\varepsilon = 0\f$ the exact
 * diode is recovered and the Newton iteration becomes a semi-smooth
 * (active-set) iteration; for \f$\varepsilon > 0\f$ the equation is smooth
 * and the valve leaks a backward flow of at most
 * \f$\varepsilon / R\f$ at \f$\Delta p = 0\f$, decaying like
 * \f$\varepsilon^2 / (R\,|\Delta p|)\f$ when the valve is firmly closed.
 *
 * Compared with \ref ValveTanh this formulation
 *
 * * has only one open resistance instead of an open/closed resistance pair,
 *   so there is no \f$R_{max}/R_{min}\f$ ratio to condition the system,
 * * has a Jacobian whose entries are bounded by \f$1/2\f$ and \f$R\f$ for
 *   every state, instead of growing with
 *   \f$k\,Q\,(R_{max}-R_{min})\f$ at the transition,
 * * converges to a well-defined limit (the ideal diode) as
 *   \f$\varepsilon \to 0\f$, so tightening the regularization refines the
 *   solution instead of changing it qualitatively.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [D, l=$R$, *-*] (3,0)
 * node[anchor=south]{$P_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * RQ_{in} - \frac{1}{2}\left(P_{in}-P_{out} +
 * \sqrt{(P_{in}-P_{out})^2 + 4\varepsilon^2}\right) = 0
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
 * -1/2 & R &  1/2 & 0\\
 * 0    & 1 &  0   & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * -\frac{1}{2}\sqrt{\Delta p^2 + 4\varepsilon^2} \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \left(\frac{\partial\mathbf{c}}{\partial\mathbf{y}}\right)^{e} =
 * \left[\begin{array}{cccc}
 * (1/2 - s) & 0 &  (s - 1/2) & 0\\
 * 0         & 0 &  0         & 0
 * \end{array}\right] \f]
 * with
 * \f[
 * s = \frac{1}{2}\left(1 +
 * \frac{\Delta p}{\sqrt{\Delta p^2 + 4\varepsilon^2}}\right).
 * \f]
 *
 * The valve status is a post-processing variable, see
 * DOFHandler::register_output_variable: it is reported but has no equation
 * and is computed from the converged state in \ref update_outputs.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` R: Resistance of the open valve (must be positive)
 * * `1` eps: Smoothing pressure, optional, defaults to 0 (exact diode)
 * * `2` upstream_block: Name of block connected upstream
 * * `3` downstream_block: Name of block connected downstream
 *
 * ### Usage in json configuration file
 *
 *     "valves": [
 *         {
 *             "type": "ValveDiode",
 *             "name": "valve",
 *             "params": {
 *                 "R": 0.005,
 *                 "eps": 0.0,
 *                 "upstream_block": "upstream_vessel",
 *                 "downstream_block": "downstream_vessel"
 *             }
 *         }
 *     ]
 *
 * ### Internal variables
 *
 * This block has no internal variables. It reports one post-processing
 * variable, `valve_status`, which has no equation and is computed from the
 * converged state once the time step is accepted.
 *
 */
class ValveDiode : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    R = 0,
    EPS = 1,
  };

  /**
   * @brief Construct a new ValveDiode object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ValveDiode(int id, Model* model)
      : Block(id, model, BlockType::valve_diode, BlockClass::valve,
              {{"R", InputParameter()},
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
   */
  void update_outputs(Eigen::Matrix<double, Eigen::Dynamic, 1>& y) override;

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{5, 0, 2};
};

#endif  // SVZERODSOLVER_MODEL_VALVEDIODE_HPP_
