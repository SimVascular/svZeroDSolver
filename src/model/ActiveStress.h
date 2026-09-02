// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file ActiveStress.h
 * @brief Active stress models for chamber blocks
 */

#ifndef SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_
#define SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Parameter.h"

/**
 * @brief Return type for ActiveStress::compute
 *
 * Bundles the active-stress rate coefficient and its positive part that
 * appear in the active-stress ODE of a chamber block:
 * \f[
 * \dot{\tau} + a \tau - \sigma_\text{max} a_+ = 0, \quad a_+ = \max(a, 0)
 * \f]
 */
struct ActiveStressResult {
  double act;       ///< Rate coefficient \f$a\f$ multiplying \f$\tau\f$
  double act_plus;  ///< Positive part of \f$a\f$, \f$a_+ = \max(a, 0)\f$
};

/**
 * @brief Base class for chamber active stress models
 *
 * Subclasses compute the active-stress rate coefficient and its positive
 * part from the activation value \f$f(t) \in [0, 1]\f$ produced by a
 * chamber's \ref ActivationFunction.
 */
class ActiveStress {
 public:
  /**
   * @brief Properties of the input parameters for this active stress model
   * [(name, InputParameter), ...]
   */
  const std::vector<std::pair<std::string, InputParameter>>
      input_param_properties;

  /**
   * @brief Construct an ActiveStress model
   *
   * @param props Parameter name/spec pairs for this active stress type
   */
  explicit ActiveStress(
      const std::vector<std::pair<std::string, InputParameter>>& props);

  virtual ~ActiveStress() = default;

  /**
   * @brief Compute active stress contributions at the current state
   *
   * @param f Activation value from the chamber's ActivationFunction,
   * \f$f \in [0, 1]\f$
   * @param alpha_max Maximum activation parameter \f$\alpha_\text{max}\f$
   * @param alpha_min Minimum activation parameter \f$\alpha_\text{min}\f$
   * @return ActiveStressResult
   */
  virtual ActiveStressResult compute(double f, double alpha_max,
                                     double alpha_min) const = 0;

  /**
   * @brief Set a scalar parameter value by name
   *
   * @param name Parameter name
   * @param value Parameter value
   */
  void set_param(const std::string& name, double value);

 protected:
  /**
   * @brief Map of parameter names to their values
   */
  std::map<std::string, double> params_;
};

/**
 * @brief Elastance-type active stress model
 *
 * Implements the active-stress rate coefficient used by \ref ChamberSphere:
 * \f[
 * a = \alpha_\text{max} f + \alpha_\text{min} (1 - f)
 * \f]
 */
class ElastanceActiveStress : public ActiveStress {
 public:
  ElastanceActiveStress() : ActiveStress({}) {}

  ActiveStressResult compute(double f, double alpha_max,
                             double alpha_min) const override;
};

#endif  // SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_
