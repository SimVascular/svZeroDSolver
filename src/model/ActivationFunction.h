// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file ActivationFunction.h
 * @brief Activation function classes for cardiac chamber models
 */

#ifndef SVZERODSOLVER_MODEL_ACTIVATIONFUNCTION_HPP_
#define SVZERODSOLVER_MODEL_ACTIVATIONFUNCTION_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Parameter.h"

/**
 * @brief Base class for activation functions
 *
 * Activation functions compute the activation value (between 0 and 1) at a
 * given time point within a cardiac cycle. These are used to modulate
 * chamber elastance over time.
 */
class ActivationFunction {
 public:
  /**
   * @brief Properties of the input parameters for this activation function
   * [(name, InputParameter), ...]
   */
  const std::vector<std::pair<std::string, InputParameter>>
      input_param_properties;

  /**
   * @brief Construct activation function
   *
   * @param cardiac_period Cardiac cycle period
   * @param input_param_properties Properties of the input parameters
   * [(name, InputParameter), ...] for this activation function
   */
  ActivationFunction(double cardiac_period,
                     const std::vector<std::pair<std::string, InputParameter>>&
                         input_param_properties);

  /**
   * @brief Virtual destructor
   */
  virtual ~ActivationFunction() = default;

  /**
   * @brief Compute activation value at given time
   *
   * @param time Current time
   * @return Activation value between 0 and 1
   */
  virtual double compute(double time) = 0;

  /**
   * @brief Create a default activation function from activation function type
   *
   * @param type_str One of: "half_cosine", "piecewise_cosine", "two_hill"
   * @param cardiac_period Cardiac cycle period
   * @return Unique pointer to the created activation function
   */
  static std::unique_ptr<ActivationFunction> create_default(
      const std::string& type_str, double cardiac_period);

  /**
   * @brief Set a scalar parameter value by name.
   *
   * Calling function must validate the parameter name and value
   *
   * @param name Parameter name
   * @param value Parameter value
   */
  void set_param(const std::string& name, double value);

  /**
   * @brief Called after all parameters are set (e.g. by loader).
   *
   * Default no-op. TwoHillActivation overrides to recompute normalization.
   */
  virtual void finalize() {}

 protected:
  /**
   * @brief Time duration of one cardiac cycle
   */
  double cardiac_period_;

  /**
   * @brief Map of parameter names to their values
   */
  std::map<std::string, double> params_;
};

/**
 * @brief Half cosine activation function
 *
 * This implements the activation function used in the original
 * ChamberElastanceInductor. The activation follows a half cosine wave
 * during the contraction period.
 *
 * \f[
 * A(t) = \begin{cases}
 * -\frac{1}{2}\cos(2\pi t_{contract}/t_{twitch}) + \frac{1}{2}, & \text{if }
 * t_{contract} \le t_{twitch} \\ 0, & \text{otherwise}
 * \end{cases}
 * \f]
 *
 * where \f$t_{contract} = \max(0, t_{in\_cycle} - t_{active})\f$
 */
class HalfCosineActivation : public ActivationFunction {
 public:
  /**
   * @brief Construct with default parameter values (loader fills via
   * set_param).
   *
   * @param cardiac_period Cardiac cycle period
   */
  explicit HalfCosineActivation(double cardiac_period)
      : ActivationFunction(cardiac_period, {{"t_active", InputParameter()},
                                            {"t_twitch", InputParameter()}}) {}

  double compute(double time) override;
};

/**
 * @brief Piecewise cosine activation function
 *
 * This implements the activation function from the LinearElastanceChamber
 * (Regazzoni chamber model). The activation consists of separate contraction
 * and relaxation phases, each following a cosine curve.
 *
 * \f[
 * \phi(t, t_C, t_R, T_C, T_R) = \begin{cases}
 * \frac{1}{2}\left[1 - \cos\left(\frac{\pi}{T_C} \bmod(t - t_C,
 * T_{HB})\right)\right],
 *   & \text{if } 0 \le \bmod(t - t_C, T_{HB}) < T_C \\
 * \frac{1}{2}\left[1 + \cos\left(\frac{\pi}{T_R} \bmod(t - t_R,
 * T_{HB})\right)\right],
 *   & \text{if } 0 \le \bmod(t - t_R, T_{HB}) < T_R \\
 * 0, & \text{otherwise}
 * \end{cases}
 * \f]
 */
class PiecewiseCosineActivation : public ActivationFunction {
 public:
  /**
   * @brief Construct with default parameter values (loader fills via
   * set_param).
   *
   * @param cardiac_period Cardiac cycle period
   */
  explicit PiecewiseCosineActivation(double cardiac_period)
      : ActivationFunction(cardiac_period,
                           {{"contract_start", InputParameter()},
                            {"relax_start", InputParameter()},
                            {"contract_duration", InputParameter()},
                            {"relax_duration", InputParameter()}}) {}

  double compute(double time) override;
};

/**
 * @brief Two hill activation function
 *
 * This implements the two-hill activation function which provides more
 * flexible and physiologically realistic waveforms. See
 * https://link.springer.com/article/10.1007/s10439-022-03047-3
 *
 * The activation is computed as:
 * \f[
 * A(t) = C \cdot \frac{g_1(t)}{1 + g_1(t)} \cdot \frac{1}{1 + g_2(t)}
 * \f]
 *
 * where:
 * \f[
 * g_1(t) = \left(\frac{t_{shifted}}{\tau_1}\right)^{m_1}, \quad
 * g_2(t) = \left(\frac{t_{shifted}}{\tau_2}\right)^{m_2}
 * \f]
 *
 * and \f$t_{shifted} = (t - t_{shift}) \bmod T_{cardiac}\f$, and \f$C\f$ is a
 * normalization constant to ensure max activation is 1.
 */
class TwoHillActivation : public ActivationFunction {
 public:
  /**
   * @brief Construct with default parameter values (loader fills via
   * set_param).
   *
   * @param cardiac_period Cardiac cycle period
   */
  explicit TwoHillActivation(double cardiac_period)
      : ActivationFunction(cardiac_period, {{"t_shift", InputParameter()},
                                            {"tau_1", InputParameter()},
                                            {"tau_2", InputParameter()},
                                            {"m1", InputParameter()},
                                            {"m2", InputParameter()}}),
        normalization_factor_(1.0),
        normalization_initialized_(false) {}

  double compute(double time) override;

  void finalize() override;

 private:
  void calculate_normalization_factor();

  double normalization_factor_;
  bool normalization_initialized_;
};

#endif  // SVZERODSOLVER_MODEL_ACTIVATIONFUNCTION_HPP_
