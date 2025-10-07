// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ClosedLoopCoronaryRightBC.h
 * @brief Right side of ClosedLoopCoronaryBC
 */
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYRIGHTBC_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYRIGHTBC_HPP_

#include "ClosedLoopCoronaryBC.h"

/**
 * @brief Right side of closed loop coronary boundary condition
 * ClosedLoopCoronaryBC.
 */
class ClosedLoopCoronaryRightBC : public ClosedLoopCoronaryBC {
 public:
  /**
   * @brief Construct a new ClosedLoopCoronaryRightBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ClosedLoopCoronaryRightBC(int id, Model *model)
      : ClosedLoopCoronaryBC(id, model,
                             BlockType::closed_loop_coronary_right_bc) {}

  /**
   * @brief Setup parameters that depend on the model
   *
   */
  void setup_model_dependent_params();
};

#endif  // SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYRIGHTBC_HPP_
