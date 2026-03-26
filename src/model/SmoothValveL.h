// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#ifndef SVZERODSOLVER_MODEL_SMOOTHVALVEL_HPP_
#define SVZERODSOLVER_MODEL_SMOOTHVALVEL_HPP_

#include <cmath>
#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Smooth valve with inductance matching the monolithic outflow equation.
 *
 * Combines L and R in one equation: P_in - P_out - R_eff*Q - L*dQ = 0
 * where R_eff = Rmin*s + Rmax*(1-s), s = 0.5*(1+tanh(k*ΔP)).
 *
 * Open (s=1): L*dQ + Rmin*Q = ΔP — correct valve resistance with inertia
 * Closed (s=0): L*dQ + Rmax*Q = ΔP — high resistance blocks backflow
 *
 * The R_eff*Q product is in C (not F) with dC/dQ = R_eff. The valve sigmoid
 * is treated as constant during linearization (no dC/dP), matching the
 * monolithic's treatment of valve coefficients.
 */
class SmoothValveL : public Block {
 public:
  enum ParamId { RMIN = 0, RMAX = 1, INDUCTANCE = 2, STEEPNESS = 3 };

  SmoothValveL(int id, Model* model)
      : Block(id, model, BlockType::smooth_valve_l, BlockClass::valve,
              {{"Rmin", InputParameter()},
               {"Rmax", InputParameter()},
               {"Impedance", InputParameter()},
               {"Steepness", InputParameter()},
               {"upstream_block", InputParameter(false, false, false)},
               {"downstream_block", InputParameter(false, false, false)}}) {}

  void setup_dofs(DOFHandler& dofhandler);
  void update_constant(SparseSystem& system, std::vector<double>& parameters);
  void update_solution(SparseSystem& system, std::vector<double>& parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy);

  TripletsContributions num_triplets{4, 1, 1};
};

#endif
