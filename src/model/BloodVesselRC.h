// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#ifndef SVZERODSOLVER_MODEL_BLOODVESSELRC_HPP_
#define SVZERODSOLVER_MODEL_BLOODVESSELRC_HPP_

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Flow-through RC windkessel for pulmonary circulation.
 *
 * Models a resistance-capacitance windkessel where the inlet and outlet
 * flows are equal (Q_in = Q_out). The capacitor stores pressure but does
 * not store blood volume. This matches the pulmonary model in
 * \cite sankaran2012patient (equation X6' = X4 - X7).
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
  /**
   * @brief Construct a new BloodVesselRC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  BloodVesselRC(int id, Model* model)
      : Block(id, model, BlockType::blood_vessel_rc, BlockClass::vessel,
              {{"Rpd", InputParameter()}, {"Cp", InputParameter()}}) {}

  /// @brief Local IDs of the parameters
  enum ParamId { RPD = 0, CP = 1 };

  void setup_dofs(DOFHandler& dofhandler);
  void update_constant(SparseSystem& system, std::vector<double>& parameters);

  /// @brief Number of triplets of element
  TripletsContributions num_triplets{5, 1, 0};
};

#endif  // SVZERODSOLVER_MODEL_BLOODVESSELRC_HPP_
