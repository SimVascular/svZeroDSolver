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
 * Inherits Emax, Emin, Vrd, Vrs, Impedance from ChamberElastanceInductor.
 * Additional parameters:
 *
 * * `5` Kxp: Passive pressure scaling
 * * `6` Kxv: Passive volume scaling
 * * `7` Vaso: Passive resting volume
 *
 */
class ChamberElastanceInductorExponential : public ChamberElastanceInductor {
 public:
  ChamberElastanceInductorExponential(int id, Model* model)
      : ChamberElastanceInductor(
            id, model,
            BlockType::chamber_elastance_inductor_exponential,
            {{"Emax", InputParameter()},
             {"Emin", InputParameter()},
             {"Vrd", InputParameter()},
             {"Vrs", InputParameter()},
             {"Impedance", InputParameter()},
             {"Kxp", InputParameter()},
             {"Kxv", InputParameter()},
             {"Vaso", InputParameter()}}) {}

  enum ExponentialParamId {
    KXP = 5,
    KXV = 6,
    VASO = 7,
  };

  void update_time(SparseSystem& system, std::vector<double>& parameters);
  void update_solution(SparseSystem& system, std::vector<double>& parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy);

  TripletsContributions num_triplets{6, 2, 1};

 protected:
  void get_elastance_values(std::vector<double>& parameters) override;
};

#endif  // SVZERODSOLVER_MODEL_CHAMBERELASTANCEINDUCTOREXPONENTIAL_HPP_
