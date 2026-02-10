// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file ActivationFunction.h
 * @brief Activation function classes for cardiac chamber models
 */

#ifndef SVZERODSOLVER_MODEL_ACTIVATIONFUNCTION_HPP_
#define SVZERODSOLVER_MODEL_ACTIVATIONFUNCTION_HPP_

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @brief Enum for activation function types
 */
enum class ActivationType {
  HALF_COSINE = 0,
  PIECEWISE_COSINE = 1,
  TWO_HILL = 2
};

/**
 * @brief Parameters for creating activation functions
 */
struct ActivationFunctionParams {
  ActivationType type;
  double cardiac_period;
  
  // Half cosine parameters
  double t_active = 0.0;
  double t_twitch = 0.0;
  
  // Piecewise cosine parameters
  double contract_start = 0.0;
  double relax_start = 0.0;
  double contract_duration = 0.0;
  double relax_duration = 0.0;
  
  // Two hill parameters
  double t_shift = 0.0;
  double tau_1 = 0.0;
  double tau_2 = 0.0;
  double m1 = 0.0;
  double m2 = 0.0;
};

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
   * @brief Virtual destructor
   */
  virtual ~ActivationFunction() = default;

  /**
   * @brief Compute activation value at given time
   * 
   * @param time Current time
   * @param cardiac_period Period of cardiac cycle
   * @return Activation value between 0 and 1
   */
  virtual double compute(double time, double cardiac_period) = 0;

  /**
   * @brief Get the type of activation function
   * 
   * @return Type of activation function
   */
  virtual ActivationType get_type() const = 0;
  
  /**
   * @brief Factory method to create activation function from parameters
   * 
   * @param params Parameters for creating the activation function
   * @return Unique pointer to the created activation function
   */
  static std::unique_ptr<ActivationFunction> create(
      const ActivationFunctionParams& params);
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
 * -\frac{1}{2}\cos(2\pi t_{contract}/t_{twitch}) + \frac{1}{2}, & \text{if } t_{contract} \le t_{twitch} \\
 * 0, & \text{otherwise}
 * \end{cases}
 * \f]
 * 
 * where \f$t_{contract} = \max(0, t_{in\_cycle} - t_{active})\f$
 */
class HalfCosineActivation : public ActivationFunction {
 public:
  /**
   * @brief Construct a new Half Cosine Activation object
   * 
   * @param t_active Time when activation begins within cardiac cycle
   * @param t_twitch Duration of the contraction twitch
   */
  HalfCosineActivation(double t_active, double t_twitch)
      : t_active_(t_active), t_twitch_(t_twitch) {}

  /**
   * @brief Compute activation value
   * 
   * @param time Current time
   * @param cardiac_period Period of cardiac cycle
   * @return Activation value between 0 and 1
   */
  double compute(double time, double cardiac_period) override {
    double t_in_cycle = std::fmod(time, cardiac_period);
    
    double t_contract = 0.0;
    if (t_in_cycle >= t_active_) {
      t_contract = t_in_cycle - t_active_;
    }

    double act = 0.0;
    if (t_contract <= t_twitch_) {
      act = -0.5 * std::cos(2.0 * M_PI * t_contract / t_twitch_) + 0.5;
    }

    return act;
  }

  /**
   * @brief Get the type of activation function
   * 
   * @return ActivationType::HALF_COSINE
   */
  ActivationType get_type() const override {
    return ActivationType::HALF_COSINE;
  }

 private:
  double t_active_;  ///< Time when activation begins
  double t_twitch_;  ///< Duration of contraction twitch
};

/**
 * @brief Piecewise cosine activation function
 * 
 * This implements the activation function from the PiecewiseCosineChamber
 * (Regazzoni chamber model). The activation consists of separate contraction
 * and relaxation phases, each following a cosine curve.
 * 
 * \f[
 * \phi(t, t_C, t_R, T_C, T_R) = \begin{cases}
 * \frac{1}{2}\left[1 - \cos\left(\frac{\pi}{T_C} \bmod(t - t_C, T_{HB})\right)\right], 
 *   & \text{if } 0 \le \bmod(t - t_C, T_{HB}) < T_C \\
 * \frac{1}{2}\left[1 + \cos\left(\frac{\pi}{T_R} \bmod(t - t_R, T_{HB})\right)\right], 
 *   & \text{if } 0 \le \bmod(t - t_R, T_{HB}) < T_R \\
 * 0, & \text{otherwise}
 * \end{cases}
 * \f]
 */
class PiecewiseCosineActivation : public ActivationFunction {
 public:
  /**
   * @brief Construct a new Piecewise Cosine Activation object
   * 
   * @param contract_start Time when contraction starts
   * @param relax_start Time when relaxation starts
   * @param contract_duration Duration of contraction phase
   * @param relax_duration Duration of relaxation phase
   */
  PiecewiseCosineActivation(double contract_start, double relax_start,
                            double contract_duration, double relax_duration)
      : contract_start_(contract_start),
        relax_start_(relax_start),
        contract_duration_(contract_duration),
        relax_duration_(relax_duration) {}

  /**
   * @brief Compute activation value
   * 
   * @param time Current time
   * @param cardiac_period Period of cardiac cycle
   * @return Activation value between 0 and 1
   */
  double compute(double time, double cardiac_period) override {
    double phi = 0.0;

    double piecewise_condition = std::fmod(time - contract_start_, cardiac_period);

    if (0.0 <= piecewise_condition && piecewise_condition < contract_duration_) {
      phi = 0.5 * (1.0 - std::cos((M_PI * piecewise_condition) / contract_duration_));
    } else {
      piecewise_condition = std::fmod(time - relax_start_, cardiac_period);
      if (0.0 <= piecewise_condition && piecewise_condition < relax_duration_) {
        phi = 0.5 * (1.0 + std::cos((M_PI * piecewise_condition) / relax_duration_));
      }
    }

    return phi;
  }

  /**
   * @brief Get the type of activation function
   * 
   * @return ActivationType::PIECEWISE_COSINE
   */
  ActivationType get_type() const override {
    return ActivationType::PIECEWISE_COSINE;
  }

 private:
  double contract_start_;     ///< Start time of contraction
  double relax_start_;        ///< Start time of relaxation
  double contract_duration_;  ///< Duration of contraction
  double relax_duration_;     ///< Duration of relaxation
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
   * @brief Construct a new Two Hill Activation object
   * 
   * @param t_shift Time shift parameter
   * @param tau_1 Time constant for first hill
   * @param tau_2 Time constant for second hill
   * @param m1 Exponent for first hill
   * @param m2 Exponent for second hill
   * @param cardiac_period Cardiac cycle period (needed for normalization)
   */
  TwoHillActivation(double t_shift, double tau_1, double tau_2, double m1,
                    double m2, double cardiac_period)
      : t_shift_(t_shift),
        tau_1_(tau_1),
        tau_2_(tau_2),
        m1_(m1),
        m2_(m2),
        normalization_factor_(1.0),
        normalization_initialized_(false) {
    initialize_normalization(cardiac_period);
  }

  /**
   * @brief Compute activation value
   * 
   * @param time Current time
   * @param cardiac_period Period of cardiac cycle
   * @return Activation value between 0 and 1
   */
  double compute(double time, double cardiac_period) override {
    if (!normalization_initialized_) {
      throw std::runtime_error("TwoHillActivation: Normalization not initialized");
    }

    double t_in_cycle = std::fmod(time, cardiac_period);
    
    // Compute shifted time (handle negative modulo)
    double t_shifted = std::fmod(t_in_cycle - t_shift_, cardiac_period);
    t_shifted = (t_shifted >= 0.0) ? t_shifted : t_shifted + cardiac_period;

    // Compute hill functions
    double g1 = std::pow(t_shifted / tau_1_, m1_);
    double g2 = std::pow(t_shifted / tau_2_, m2_);

    // Compute activation with normalization
    double act = normalization_factor_ * (g1 / (1.0 + g1)) * (1.0 / (1.0 + g2));

    return act;
  }

  /**
   * @brief Get the type of activation function
   * 
   * @return ActivationType::TWO_HILL
   */
  ActivationType get_type() const override {
    return ActivationType::TWO_HILL;
  }

 private:
  /**
   * @brief Initialize normalization factor
   * 
   * Computes the maximum value of the two-hill function over one cardiac
   * cycle to normalize the activation to [0, 1].
   * 
   * @param cardiac_period Period of cardiac cycle
   */
  void initialize_normalization(double cardiac_period) {
    // Time step for numerical integration when finding maximum activation
    // Value chosen to balance accuracy and computational cost
    constexpr double NORMALIZATION_DT = 1e-5;
    
    double max_value = 0.0;

    for (double t_temp = 0.0; t_temp < cardiac_period; t_temp += NORMALIZATION_DT) {
      double g1 = std::pow(t_temp / tau_1_, m1_);
      double g2 = std::pow(t_temp / tau_2_, m2_);
      double two_hill_val = (g1 / (1.0 + g1)) * (1.0 / (1.0 + g2));

      max_value = std::max(max_value, two_hill_val);
    }

    normalization_factor_ = 1.0 / max_value;
    normalization_initialized_ = true;
  }

  double t_shift_;                    ///< Time shift parameter
  double tau_1_;                      ///< Time constant for first hill
  double tau_2_;                      ///< Time constant for second hill
  double m1_;                         ///< Exponent for first hill
  double m2_;                         ///< Exponent for second hill
  double normalization_factor_;       ///< Normalization constant
  bool normalization_initialized_;    ///< Flag for normalization
};

/**
 * @brief Factory method implementation for creating activation functions
 */
inline std::unique_ptr<ActivationFunction> ActivationFunction::create(
    const ActivationFunctionParams& params) {
  switch (params.type) {
    case ActivationType::HALF_COSINE:
      return std::make_unique<HalfCosineActivation>(
          params.t_active, params.t_twitch);
    
    case ActivationType::PIECEWISE_COSINE:
      return std::make_unique<PiecewiseCosineActivation>(
          params.contract_start, params.relax_start,
          params.contract_duration, params.relax_duration);
    
    case ActivationType::TWO_HILL:
      return std::make_unique<TwoHillActivation>(
          params.t_shift, params.tau_1, params.tau_2,
          params.m1, params.m2, params.cardiac_period);
    
    default:
      throw std::runtime_error(
          "ActivationFunction::create: Invalid activation type " +
          std::to_string(static_cast<int>(params.type)));
  }
}

#endif  // SVZERODSOLVER_MODEL_ACTIVATIONFUNCTION_HPP_
