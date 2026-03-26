// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#ifndef SVZERODSOLVER_MODEL_BLOODVESSELRC_HPP_
#define SVZERODSOLVER_MODEL_BLOODVESSELRC_HPP_

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Flow-through pulmonary windkessel matching the pulmonary equation
 * from ClosedLoopHeartPulmonary.
 *
 * Models a resistance-capacitance windkessel where the inlet and outlet
 * flows are equal (Q_in = Q_out). The capacitor stores pressure but does
 * not store blood volume. This matches the monolithic block where Q_RV
 * flows through the pulmonary directly to the LA volume equation.
 *
 * ### Governing equations
 *
 * \f[
 * C_p \, \dot{P}_{in} + \frac{P_{in} - P_{out}}{R_{pd}} - Q_{in} = 0
 * \f]
 * \f[
 * Q_{in} - Q_{out} = 0
 * \f]
 *
 * ### Parameters
 *
 * * `Rpd` — Pulmonary resistance
 * * `Cp` — Pulmonary capacitance
 */
class BloodVesselRC : public Block {
 public:
  BloodVesselRC(int id, Model* model)
      : Block(id, model, BlockType::blood_vessel_rc, BlockClass::vessel,
              {{"Rpd", InputParameter()}, {"Cp", InputParameter()}}) {}

  enum ParamId { RPD = 0, CP = 1 };

  void setup_dofs(DOFHandler& dofhandler);
  void update_constant(SparseSystem& system, std::vector<double>& parameters);

  TripletsContributions num_triplets{5, 1, 0};
};

#endif  // SVZERODSOLVER_MODEL_BLOODVESSELRC_HPP_
