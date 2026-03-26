// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#ifndef SVZERODSOLVER_MODEL_IDEALVALVE_HPP_
#define SVZERODSOLVER_MODEL_IDEALVALVE_HPP_

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Ideal diode valve matching the ClosedLoopHeartPulmonary valve logic.
 *
 * Combines the equation structure of PiecewiseValve (piecewise resistance for
 * stable Newton convergence) with the valve logic and post_solve enforcement
 * from the monolithic heart block.
 *
 * Differences from PiecewiseValve:
 *   - Checks both pressure AND flow direction: valve closes when
 *     P_in <= P_out **and** Q_in <= 0 (monolithic condition).
 *   - post_solve zeros Q_in and Q_out when the valve is closed, eliminating
 *     residual backflow through the finite Rmax.
 *
 * Differences from ValveTanh:
 *   - No smooth tanh transition (piecewise switch, no Jacobian contribution
 *     from the valve function itself).
 *   - Checks flow direction in addition to pressure.
 *   - post_solve zeroing.
 *
 * ### Governing equations
 *
 * \f[
 * P_{in} - P_{out} - R(v) \, Q_{in} = 0
 * \f]
 * \f[
 * Q_{in} - Q_{out} = 0
 * \f]
 *
 * where \f$ R = R_{min} \f$ when open, \f$ R = R_{max} \f$ when closed.
 * After solving, if closed: \f$ Q_{in} = Q_{out} = 0 \f$.
 *
 * ### Parameters
 *
 * * `Rmin` — Open-valve resistance
 * * `Rmax` — Closed-valve resistance (large value for Newton stability)
 */
class IdealValve : public Block {
 public:
  enum ParamId { RMIN = 0, RMAX = 1, IMPEDANCE = 2 };

  IdealValve(int id, Model* model)
      : Block(id, model, BlockType::ideal_valve, BlockClass::valve,
              {{"Rmin", InputParameter()},
               {"Rmax", InputParameter()},
               {"Impedance", InputParameter(true)},
               {"upstream_block", InputParameter(false, false, false)},
               {"downstream_block", InputParameter(false, false, false)}}) {}

  void setup_dofs(DOFHandler& dofhandler);
  void update_constant(SparseSystem& system, std::vector<double>& parameters);
  void update_solution(SparseSystem& system, std::vector<double>& parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy);
  TripletsContributions num_triplets{5, 1, 3};

 private:
  double valve_ = 1.0;
};

#endif  // SVZERODSOLVER_MODEL_IDEALVALVE_HPP_
