// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "Parameter.h"

// Only scalar expressions of t are needed; disabling the rest of ExprTk
// roughly halves the compile time and size of this file's object.
#define exprtk_disable_string_capabilities
#define exprtk_disable_rtl_io
#define exprtk_disable_rtl_io_file
#define exprtk_disable_rtl_vecops
#define exprtk_disable_enhanced_features
#define exprtk_disable_superscalar_unroll
#include <exprtk.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

/// Number of midpoint samples used to average an expression over one cardiac
/// cycle when converting it to a steady value.
constexpr int kNumSteadyMeanSamples = 1000;

/// PIMPL implementation holding the exprtk expression state
struct Parameter::ExprtkImpl {
  exprtk::symbol_table<double>
      symbol_table;  ///< Symbol table mapping 't' to current time
  exprtk::expression<double> expression;  ///< Compiled expression object
  double* time_value = nullptr;  ///< Pointer into the symbol table for 't'
};

Parameter::~Parameter() = default;
Parameter::Parameter(Parameter&&) noexcept = default;
Parameter& Parameter::operator=(Parameter&&) noexcept = default;

namespace {
// Checks that (xs, ys) form a valid piecewise-linear break point table for
// linear_interpolate/linear_interpolate_slope. xs must be strictly
// increasing; a near-vertical segment (two xs values an epsilon apart) can
// be used instead of a literal duplicate to approximate a step change.
void validate_breakpoint_table(const std::vector<double>& xs,
                               const std::vector<double>& ys) {
  if (xs.size() < 2) {
    throw std::runtime_error(
        "linear_interpolate: break point table must have at least 2 "
        "entries (got " +
        std::to_string(xs.size()) + ")");
  }
  if (xs.size() != ys.size()) {
    throw std::runtime_error(
        "linear_interpolate: xs and ys must have the same length (got " +
        std::to_string(xs.size()) + " and " + std::to_string(ys.size()) + ")");
  }
  for (size_t i = 0; i + 1 < xs.size(); i++) {
    if (xs[i] >= xs[i + 1]) {
      throw std::runtime_error(
          "linear_interpolate: xs must be strictly increasing (xs[" +
          std::to_string(i) + "] = " + std::to_string(xs[i]) +
          " is not less than xs[" + std::to_string(i + 1) +
          "] = " + std::to_string(xs[i + 1]) + ")");
    }
  }
}
}  // namespace

Parameter::Parameter(int id, double value) {
  this->id = id;
  update(value);
}

Parameter::Parameter(int id, const std::vector<double>& times,
                     const std::vector<double>& values, bool periodic) {
  this->id = id;
  this->is_periodic = periodic;
  update(times, values);
}

Parameter::Parameter(int id, const std::string& expression_string) {
  this->id = id;
  update(expression_string);
}

void Parameter::update(double update_value) {
  is_function = false;
  exprtk_impl_.reset();
  is_constant = true;
  is_periodic = true;
  value = update_value;
}

void Parameter::update(const std::vector<double>& update_times,
                       const std::vector<double>& update_values) {
  is_function = false;
  exprtk_impl_.reset();
  this->size = update_values.size();

  if (size == 1) {
    value = update_values[0];
    is_constant = true;
  } else {
    times = update_times;
    values = update_values;
    cycle_period = update_times.back() - update_times[0];
    is_constant = false;
  }
}

void Parameter::update(const std::string& update_string) {
  auto impl = std::make_unique<ExprtkImpl>();
  impl->symbol_table.create_variable("t");
  impl->symbol_table.add_constants();  // pi, epsilon, inf
  impl->time_value = &impl->symbol_table.get_variable("t")->ref();
  impl->expression.register_symbol_table(impl->symbol_table);

  // Reject loops, return statements, assignments and variable definitions:
  // a boundary condition only needs a value of t, and a loop could hang the
  // solver while assignments could silently modify t.
  using settings_t = exprtk::parser<double>::settings_t;
  settings_t settings;
  settings.disable_control_structure(settings_t::e_ctrl_for_loop)
      .disable_control_structure(settings_t::e_ctrl_while_loop)
      .disable_control_structure(settings_t::e_ctrl_repeat_loop)
      .disable_control_structure(settings_t::e_ctrl_return)
      .disable_all_assignment_ops()
      .disable_local_vardef();

  exprtk::parser<double> parser(settings);
  if (!parser.compile(update_string, impl->expression)) {
    std::ostringstream msg;
    msg << "could not compile expression \"" << update_string << "\"";
    for (std::size_t i = 0; i < parser.error_count(); ++i) {
      const auto error = parser.get_error(i);
      msg << "; ";
      // Errors that are not tied to a token report position npos
      if (error.token.position != std::string::npos) {
        msg << "at position " << error.token.position << ": ";
      }
      msg << error.diagnostic;
    }
    throw std::runtime_error(msg.str());
  }

  // Only modify this parameter once the expression is known to be valid.
  exprtk_impl_ = std::move(impl);
  is_function = true;
  is_constant = false;
  is_periodic = false;
}

double Parameter::get(double time) {
  // Return the constant value if parameter is constant
  if (is_constant) {
    return value;
  }

  if (is_function) {
    *exprtk_impl_->time_value = time;
    return exprtk_impl_->expression.value();
  }

  // Determine the time within this->times (necessary to extrapolate)
  double rtime;

  if (is_periodic == true) {
    rtime = fmod(time, cycle_period);
  } else {
    // this->times is not periodic when running with external solver
    rtime = time;
  }

  // Determine the lower and upper element for interpolation
  auto i = lower_bound(times.begin(), times.end(), rtime);
  int k = i - times.begin();

  if (i == times.end()) {
    --i;
  } else if (*i == rtime) {
    return values[k];
  }
  int m = k ? k - 1 : 1;

  // Perform linear interpolation
  // TODO: Implement periodic cubic spline
  return values[m] +
         ((values[k] - values[m]) / (times[k] - times[m])) * (rtime - times[m]);
}

void Parameter::to_steady(double cardiac_cycle_period) {
  if (is_constant) {
    return;
  }

  if (is_function) {
    // Midpoint-rule time average over one cardiac cycle
    double sum = 0.0;
    for (int i = 0; i < kNumSteadyMeanSamples; i++) {
      sum += get((i + 0.5) * cardiac_cycle_period / kNumSteadyMeanSamples);
    }
    value = sum / kNumSteadyMeanSamples;
  } else {
    value = std::accumulate(values.begin(), values.end(), 0.0) / double(size);
  }
  is_constant = true;
  steady_converted = true;
}

void Parameter::to_unsteady() {
  if (steady_converted) {
    is_constant = false;
    steady_converted = false;
  }
}

double linear_interpolate(double x, const std::vector<double>& xs,
                          const std::vector<double>& ys) {
  validate_breakpoint_table(xs, ys);

  if (x <= xs.front()) {
    return ys.front();
  }
  if (x >= xs.back()) {
    return ys.back();
  }

  auto i = lower_bound(xs.begin(), xs.end(), x);
  int k = i - xs.begin();
  if (*i == x) {
    return ys[k];
  }
  int m = k - 1;

  return ys[m] + ((ys[k] - ys[m]) / (xs[k] - xs[m])) * (x - xs[m]);
}

double linear_interpolate_slope(double x, const std::vector<double>& xs,
                                const std::vector<double>& ys) {
  validate_breakpoint_table(xs, ys);

  if (x <= xs.front() || x >= xs.back()) {
    return 0.0;
  }

  auto i = lower_bound(xs.begin(), xs.end(), x);
  int k = i - xs.begin();
  if (*i == x) {
    ++k;
  }
  int m = k - 1;

  return (ys[k] - ys[m]) / (xs[k] - xs[m]);
}
