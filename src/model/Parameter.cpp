// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "Parameter.h"

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

  exprtk::parser<double> parser;
  if (!parser.compile(update_string, impl->expression)) {
    std::ostringstream msg;
    msg << "could not compile expression \"" << update_string << "\"";
    for (std::size_t i = 0; i < parser.error_count(); ++i) {
      const auto error = parser.get_error(i);
      msg << "; at position " << error.token.position << ": "
          << error.diagnostic;
    }
    throw std::runtime_error(msg.str());
  }

  // Only modify this parameter once the expression is known to be valid.
  exprtk_impl_ = std::move(impl);
  expression_string = update_string;
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
