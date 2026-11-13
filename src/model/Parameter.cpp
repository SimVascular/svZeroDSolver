// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "Parameter.h"

#include <stdexcept>
#include <string>

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

void Parameter::update(double update_value) {
  is_constant = true;
  is_periodic = true;
  value = update_value;
}

void Parameter::update(const std::vector<double>& update_times,
                       const std::vector<double>& update_values) {
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
