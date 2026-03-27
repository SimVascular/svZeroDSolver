// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ChamberElastanceInductorExponential.h
 * @brief model::ChamberElastanceInductorExponential source file
 */
#ifndef SVZERODSOLVER_MODEL_CHAMBERELASTANCEINDUCTOREXPONENTIAL_HPP_
#define SVZERODSOLVER_MODEL_CHAMBERELASTANCEINDUCTOREXPONENTIAL_HPP_

#include "ChamberElastanceInductor.h"

/**
 * @brief Cardiac chamber with exponential passive pressure-volume relation.
 *
 * Extends ChamberElastanceInductor with an exponential passive P-V
 * relationship for modeling atrial chambers. Based on the atrial model in
 * \cite sankaran2012patient and \cite menon2023predictors.
 *
 * ### Governing equations
 *
 * The pressure-volume relation replaces the linear form with:
 *
 * \f[
 * P_{in} = A(t) \, E_{max} (V_c - V_{aso})
 *        + (1 - A(t)) \, K_{xp} \left( e^{K_{xv}(V_c - V_{aso})} - 1 \right)
 * \f]
 *
 * The inductor and volume conservation equations are inherited unchanged.
 *
 * ### Parameters
 *
 * ### Parameters
 *
 * * `0` Emax: Maximum (active) elastance
 * * `1` Impedance: Outflow inductance
 * * `2` Kxp: Passive pressure scaling
 * * `3` Kxv: Passive volume scaling
 * * `4` Vaso: Passive resting volume
 *
 */
class ChamberElastanceInductorExponential : public ChamberElastanceInductor {
 public:
  /**
   * @brief Construct a new ChamberElastanceInductorExponential object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ChamberElastanceInductorExponential(int id, Model* model)
      : ChamberElastanceInductor(
            id, model, BlockType::chamber_elastance_inductor_exponential,
            {{"Impedance", InputParameter()},
             {"Emax", InputParameter()},
             {"Kxp", InputParameter()},
             {"Kxv", InputParameter()},
             {"Vaso", InputParameter()}}) {}

  /**
   * @brief Local IDs of the parameters (Impedance=0 and Emax=1 shared with
   * base class)
   */
  enum ExponentialParamId {
    KXP = 2,
    KXV = 3,
    VASO = 4,
  };

  void update_time(SparseSystem& system, std::vector<double>& parameters);
  void update_solution(SparseSystem& system, std::vector<double>& parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy);

  /// @brief Number of triplets of element
  TripletsContributions num_triplets{6, 2, 1};

 protected:
  void get_elastance_values(std::vector<double>& parameters) override;
};

#endif  // SVZERODSOLVER_MODEL_CHAMBERELASTANCEINDUCTOREXPONENTIAL_HPP_
