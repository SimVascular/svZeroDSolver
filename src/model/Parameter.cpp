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
  time_value = nullptr;

  expression.release();
  symbol_table.clear();

  if (!symbol_table.create_variable("t")) {
    throw std::runtime_error("Error failed to create time_value in symbol_table.");
    return;
  }

  symbol_table.add_constants();

  time_value = &symbol_table.get_variable("t")->ref();
  expression.register_symbol_table(symbol_table);

  exprtk::parser<double> parser;

  if (!parser.compile(expression_string, expression)) {
    is_function = false;

    printf("Error: %s\tExpression: %s\n", parser.error().c_str(),
           expression_string.c_str());

    for (std::size_t i = 0; i < parser.error_count(); ++i) {
      typedef exprtk::parser_error::type err_t;
      const auto error = parser.get_error(i);

      printf(
          "Error: %02d  Position: %02d Type: [%14s] Msg: %s\tExpression: %s\n",
          static_cast<unsigned int>(i),
          static_cast<unsigned int>(error.token.position),
          exprtk::parser_error::to_str(error.mode).c_str(),
          error.diagnostic.c_str(), expression_string.c_str());
    }
    throw std::runtime_error("Error when compiling the function provided in 'fn'.");
    return;
  }
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
    assert(time_value != nullptr);
    *time_value = time;
    return expression.value();
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
