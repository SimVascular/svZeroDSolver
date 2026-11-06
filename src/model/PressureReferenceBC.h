// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file PressureReferenceBC.h
 * @brief model::PressureReferenceBC source file
 */
#ifndef SVZERODSOLVER_MODEL_PRESSUREREFERENCEBC_HPP_
#define SVZERODSOLVER_MODEL_PRESSUREREFERENCEBC_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief Pressure reference boundary condition.
 *
 * Applies a predefined pressure at a boundary.
 *
 * \f[
 * \begin{circuitikz}
 * \draw (1,0) node[anchor=south]{$P$} to [short, *-] (1.2,0) ;
 * \draw [-latex] (1.4,0) -- (2.2,0) node[right] {$Q$};
 * \draw (1,0) to [short, l=, *-] (1,-1)
 * node[ground]{$\hat{P}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P=\hat{P}
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{ll}P & Q\end{array}\right]^{T}
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ll}1 & 0\end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{C}^{e}=\left[\hat{P}\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Pressure
 *
 * `"P"` is given either as an array of values at the times in the array
 * `"t"`, linearly interpolated and repeated every cardiac cycle, or as a
 * string with an expression of the simulation time `t` (then `"t"` must be
 * omitted). Expressions are parsed with ExprTk
 * (https://www.partow.net/programming/exprtk/), which supports the usual
 * operators, functions (`sin`, `cos`, `exp`, `log`, ...), the constant `pi`
 * and conditionals (`t < 0.3 ? 1 : 0`); loops and assignments are disabled.
 * They are evaluated at the absolute simulation time, so a pulsatile
 * expression should be periodic with the model's cardiac cycle for the
 * simulation to reach a periodic state. That period is set by the model's
 * tabular parameters (e.g. other boundary conditions given as arrays), or
 * else by `"cardiac_period"` in `"simulation_parameters"`, or else defaults
 * to 1.0. For the steady initial condition, an expression is replaced by its
 * average over one cardiac cycle.
 *
 * ### Usage in json configuration file
 *
 *     "boundary_conditions": [
 *         {
 *             "bc_name": "OUT",
 *             "bc_type": "PRESSURE",
 *             "bc_values": {
 *                 "P": [
 *                     1000.0,
 *                     1000.0
 *                 ],
 *                 "t": [
 *                     0.0,
 *                     1.0
 *                 ]
 *             }
 *         }
 *     ]
 *
 * or, with an expression of time:
 *
 *     "boundary_conditions": [
 *         {
 *             "bc_name": "OUT",
 *             "bc_type": "PRESSURE",
 *             "bc_values": {
 *                 "P": "1000.0 + 100.0 * cos(2 * pi * t)"
 *             }
 *         }
 *     ]
 *
 * ### Internal variables
 *
 * This block has no internal variables.
 *
 */
class PressureReferenceBC : public Block {
 public:
  /**
   * @brief Construct a new PressureReferenceBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  PressureReferenceBC(int id, Model* model)
      : Block(id, model, BlockType::pressure_bc, BlockClass::boundary_condition,
              {{"t", InputParameter(false, true)},
               {"P", InputParameter(false, true)}}) {}

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
  void setup_dofs(DOFHandler& dofhandler) override;

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem& system,
                       std::vector<double>& parameters) override;

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(SparseSystem& system,
                   std::vector<double>& parameters) override;

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{1, 0, 0};
};

#endif  // SVZERODSOLVER_MODEL_PRESSUREREFERENCEBC_HPP_
