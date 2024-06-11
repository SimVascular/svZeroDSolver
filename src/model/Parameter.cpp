// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "Parameter.h"

#include <cstdio>
#include <string>

Parameter::Parameter(int id, double value) {
  this->id = id;
  update(value);
}

Parameter::Parameter(int id, const std::vector<double> &times,
                     const std::vector<double> &values, bool periodic) {
  this->id = id;
  this->is_periodic = periodic;
  update(times, values);
}

Parameter::Parameter(int id, const std::string expression_string) {
  this->id = id;
  this->expression_string = expression_string;
  update(expression_string);
}

void Parameter::update(double update_value) {
  is_constant = true;
  is_periodic = true;
  value = update_value;
}

void Parameter::update(const std::vector<double> &update_times,
                       const std::vector<double> &update_values) {
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

void Parameter::update(const std::string update_string) {
  is_function = true;
  is_constant = false;
  expression_string = update_string;
}

double Parameter::get(double time) {
  // Return the constant value if parameter is constant
  if (is_constant) {
    return value;
  }

  // Determine the time within this->times (necessary to extrapolate)
  double rtime;

  if (is_periodic == true) {
    rtime = fmod(time, cycle_period);
  } else {
    // this->times is not periodic when running with external solver
    rtime = time;
  }

  if (is_function == true) {
    // Adapted from example from Basic Design example at
    // http://www.partow.net/programming/exprtk/index.html
    double t = time;

    exprtk::symbol_table<double> symbol_table;
    symbol_table.add_variable("t", t);

    exprtk::expression<double> expression;
    expression.register_symbol_table(symbol_table);

    exprtk::parser<double> parser;

    // Compile the parser and check that the input expression is valid
    if (!parser.compile(expression_string, expression)) {
      std::runtime_error("Error when compiling the function provided in 'fn'.");
    }
    double value = expression.value();
    return value;
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

void Parameter::to_steady() {
  if (is_constant) {
    return;
  }

  value = std::accumulate(values.begin(), values.end(), 0.0) / double(size);
  is_constant = true;
  steady_converted = true;
}

void Parameter::to_unsteady() {
  if (steady_converted) {
    is_constant = false;
    steady_converted = false;
  }
}
