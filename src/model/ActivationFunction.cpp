// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ActivationFunction.h"

#include <cmath>
#include <stdexcept>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ActivationFunction::ActivationFunction(
    double cardiac_period,
    const std::vector<std::pair<std::string, InputParameter>>&
        input_param_properties)
    : cardiac_period_(cardiac_period),
      input_param_properties(input_param_properties) {
  // Initialize parameters with default values
  for (const auto& p : this->input_param_properties) {
    if (p.second.is_number) {
      double default_val = p.second.is_optional ? p.second.default_val : 0.0;
      params_[p.first] = default_val;
    }
  }
}

void ActivationFunction::set_param(const std::string& name, double value) {
  params_[name] = value;
}

std::unique_ptr<ActivationFunction> ActivationFunction::create_default(
    const std::string& type_str, double cardiac_period) {
  if (type_str == "half_cosine") {
    return std::make_unique<HalfCosineActivation>(cardiac_period);
  }
  if (type_str == "piecewise_cosine") {
    return std::make_unique<PiecewiseCosineActivation>(cardiac_period);
  }
  if (type_str == "two_hill") {
    return std::make_unique<TwoHillActivation>(cardiac_period);
  }
  throw std::runtime_error(
      "Unknown activation_function type '" + type_str +
      "'. Must be one of: half_cosine, piecewise_cosine, two_hill");
}

double HalfCosineActivation::compute(double time) {
  double t_in_cycle = std::fmod(time, cardiac_period_);
  const double t_active = params_.at("t_active");
  const double t_twitch = params_.at("t_twitch");

  double t_contract = 0.0;
  if (t_in_cycle >= t_active) {
    t_contract = t_in_cycle - t_active;
  }

  double act = 0.0;
  if (t_contract <= t_twitch) {
    act = -0.5 * std::cos(2.0 * M_PI * t_contract / t_twitch) + 0.5;
  }

  return act;
}

double PiecewiseCosineActivation::compute(double time) {
  const double contract_start = params_.at("contract_start");
  const double relax_start = params_.at("relax_start");
  const double contract_duration = params_.at("contract_duration");
  const double relax_duration = params_.at("relax_duration");

  double phi = 0.0;
  double piecewise_condition =
      std::fmod(time - contract_start, cardiac_period_);

  if (0.0 <= piecewise_condition && piecewise_condition < contract_duration) {
    phi = 0.5 *
          (1.0 - std::cos((M_PI * piecewise_condition) / contract_duration));
  } else {
    piecewise_condition = std::fmod(time - relax_start, cardiac_period_);
    if (0.0 <= piecewise_condition && piecewise_condition < relax_duration) {
      phi =
          0.5 * (1.0 + std::cos((M_PI * piecewise_condition) / relax_duration));
    }
  }

  return phi;
}

void TwoHillActivation::calculate_normalization_factor() {
  if (cardiac_period_ <= 0.0) {
    throw std::runtime_error(
        "TwoHillActivation::calculate_normalization_factor: cardiac_period "
        "must be positive (got " +
        std::to_string(cardiac_period_) + ")");
  }

  const double tau_1 = params_.at("tau_1");
  const double tau_2 = params_.at("tau_2");
  const double m1 = params_.at("m1");
  const double m2 = params_.at("m2");

  constexpr double NORMALIZATION_DT = 1e-5;
  double max_value = 0.0;

  for (double t_temp = 0.0; t_temp < cardiac_period_;
       t_temp += NORMALIZATION_DT) {
    double g1 = std::pow(t_temp / tau_1, m1);
    double g2 = std::pow(t_temp / tau_2, m2);
    double two_hill_val = (g1 / (1.0 + g1)) * (1.0 / (1.0 + g2));
    max_value = std::max(max_value, two_hill_val);
  }

  if (!(max_value > 0.0) || !std::isfinite(max_value)) {
    throw std::runtime_error(
        "TwoHillActivation::calculate_normalization_factor: max activation "
        "value must be positive and finite (got " +
        std::to_string(max_value) +
        "). Check tau_1, tau_2, m1, m2 are valid (e.g., tau_1 > 0, tau_2 > "
        "0).");
  }

  normalization_factor_ = 1.0 / max_value;
  normalization_initialized_ = true;
}

void TwoHillActivation::finalize() { calculate_normalization_factor(); }

double TwoHillActivation::compute(double time) {
  if (!normalization_initialized_) {
    throw std::runtime_error(
        "TwoHillActivation: call finalize() after setting parameters");
  }

  const double t_shift = params_.at("t_shift");
  const double tau_1 = params_.at("tau_1");
  const double tau_2 = params_.at("tau_2");
  const double m1 = params_.at("m1");
  const double m2 = params_.at("m2");

  double t_in_cycle = std::fmod(time, cardiac_period_);
  double t_shifted = std::fmod(t_in_cycle - t_shift, cardiac_period_);
  t_shifted = (t_shifted >= 0.0) ? t_shifted : t_shifted + cardiac_period_;

  double g1 = std::pow(t_shifted / tau_1, m1);
  double g2 = std::pow(t_shifted / tau_2, m2);

  return normalization_factor_ * (g1 / (1.0 + g1)) * (1.0 / (1.0 + g2));
}
