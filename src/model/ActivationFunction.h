// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file ActivationFunction.h
 * @brief Activation function classes for cardiac chamber models
 */

#ifndef SVZERODSOLVER_MODEL_ACTIVATIONFUNCTION_HPP_
#define SVZERODSOLVER_MODEL_ACTIVATIONFUNCTION_HPP_

#include <memory>
#include <string>
#include <vector>

#include "Parameter.h"
#include "ParameterHolder.h"

/**
 * @brief Type of signal an ActivationFunction produces
 *
 * Used to catch mismatched activation_function/active_stress pairings:
 * an ActiveStress model declares the signal type it needs via
 * \ref ActiveStress::required_signal_type, and ChamberSphere::setup_dofs
 * checks it against the chamber's ActivationFunction::signal_type.
 */
enum class ActivationSignalType {
  normalized,  ///< f(t) in [0, 1], used by most activation functions
  rate         ///< Signed reaction-rate signal u(t), used by \ref
               ///< PiecewiseRateActivation
};

/**
 * @brief Base class for activation functions
 *
 * Activation functions compute the activation value (between 0 and 1) at a
 * given time point within a cardiac cycle. These are used to modulate
 * chamber elastance over time.
 */
class ActivationFunction : public ParameterHolder {
 public:
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
   * @brief Compute activation value at given time
   *
   * @param time Current time
   * @return Activation value between 0 and 1
   */
  virtual double compute(double time) = 0;

  /**
   * @brief Type of signal this activation function's compute() produces
   *
   * @return \ref ActivationSignalType::normalized for all activation
   * functions with function values between [0;1] and \ref
   * ActivationSignalType::rate for \ref PiecewiseRateActivation, which returns
   * a signed reaction-rate signal.
   */
  virtual ActivationSignalType signal_type() const {
    return ActivationSignalType::normalized;
  }

  /**
   * @brief Create a default activation function from activation function type
   *
   * @param type_str One of: "half_cosine", "piecewise_cosine", "two_hill",
   * "double_tanh", "wrapping_cosine", "fourier", "piecewise_rate"
   * @param cardiac_period Cardiac cycle period
   * @return Unique pointer to the created activation function
   */
  static std::unique_ptr<ActivationFunction> create_default(
      const std::string& type_str, double cardiac_period);

 protected:
  /**
   * @brief Time duration of one cardiac cycle
   */
  double cardiac_period_;
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

/**
 * @brief Double tanh (systole/diastole sigmoid product) activation function
 *
 * This implements the original ChamberSphere activation: a smooth indicator
 * function built from the product of two tanh sigmoids, one rising at
 * systole and one falling at diastole.
 *
 * \f[
 * f(t) = S_+ (t_{in\_cycle} - t_sys) \cdot S_- (t_{in\_cycle} - t_dias), \quad
 * S_\pm (\Delta t) = \frac{1}{2} \left(1.0 \pm
 * \text{tanh}\left( \frac{\Delta t} {\gamma}
 * \right) \right)
 * \f]
 */
class DoubleTanhActivation : public ActivationFunction {
 public:
  /**
   * @brief Construct with default parameter values (loader fills via
   * set_param).
   *
   * @param cardiac_period Cardiac cycle period
   */
  explicit DoubleTanhActivation(double cardiac_period)
      : ActivationFunction(cardiac_period, {{"tsys", InputParameter()},
                                            {"tdias", InputParameter()},
                                            {"steepness", InputParameter()}}) {}

  double compute(double time) override;
};

/**
 * @brief Wrapping cosine activation function
 *
 * Reproduces the atrial activation AA(t) from \cite sankaran2012patient
 * and \cite menon2023predictors. The cosine wraps across the cardiac
 * cycle boundary, which half_cosine and piecewise_cosine cannot do.
 *
 * Parameters:
 * * `Tsa` — Atrial systole time fraction (multiplied by cardiac period)
 * * `tpwave` — P-wave divisor (cardiac period is divided by this)
 */
class WrappingCosineActivation : public ActivationFunction {
 public:
  /**
   * @brief Construct a new WrappingCosineActivation object
   *
   * @param cardiac_period Cardiac cycle period
   */
  explicit WrappingCosineActivation(double cardiac_period)
      : ActivationFunction(cardiac_period, {{"Tsa", InputParameter()},
                                            {"tpwave", InputParameter()}}) {}

  double compute(double time) override;
};

/**
 * @brief Fourier series activation function
 *
 * 25-mode Fourier series for ventricular elastance activation, based on
 * J. Tran's tuning framework (\cite menon2023predictors). The series is
 * normalized so the output ranges from 0 (diastole) to 1 (peak systole).
 * No user parameters are needed; the Fourier coefficients are hardcoded.
 *
 * Combined with ChamberElastanceInductor using
 *   Emax = Fourier_peak * E_scale,  Emin = Fourier_trough * E_scale
 * this reproduces E(t) = Fourier(t) * E_scale exactly.
 */
class FourierActivation : public ActivationFunction {
 public:
  /**
   * @brief Construct a new FourierActivation object
   *
   * @param cardiac_period Cardiac cycle period
   */
  explicit FourierActivation(double cardiac_period)
      : ActivationFunction(cardiac_period, {}),
        norm_min_(0.0),
        norm_range_(1.0),
        normalization_initialized_(false) {}

  double compute(double time) override;
  void finalize() override;

 private:
  void calculate_normalization();
  double compute_raw(double t_in_cycle) const;

  double norm_min_;
  double norm_range_;
  bool normalization_initialized_;
};

/**
 * @brief Piecewise-linear activation rate activation signal
 *
 * Reproduces an activation rate signal \f$u(t)\f$ which is a
 * piecewise-linear function of time within the cardiac cycle, given as a
 * break point table (`u_t`/`u_values`) and evaluated with \ref
 * linear_interpolate. The sign is meaningful for the strain-dependent active
 * stress model (which splits the signal into positive/negative parts).
 *
 * `u_t` must be normalized to the cardiac cycle length, i.e.
 * \f$u_t \in [0, 1]\f$, where 0 is the start of the cycle and 1 is the end.
 * This way, the break point table does not need to be recomputed if the
 * cardiac cycle length changes.
 *
 * Parameters: `u_t`, `u_values`
 */
class PiecewiseRateActivation : public ActivationFunction {
 public:
  /**
   * @brief Construct a new PiecewiseRateActivation object
   *
   * @param cardiac_period Cardiac cycle period
   */
  explicit PiecewiseRateActivation(double cardiac_period)
      : ActivationFunction(cardiac_period,
                           {{"u_t", InputParameter(false, true)},
                            {"u_values", InputParameter(false, true)}}) {}

  double compute(double time) override;

  void finalize() override;

  ActivationSignalType signal_type() const override {
    return ActivationSignalType::rate;
  }
};

#endif  // SVZERODSOLVER_MODEL_ACTIVATIONFUNCTION_HPP_
