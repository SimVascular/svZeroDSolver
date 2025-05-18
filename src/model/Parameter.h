// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file Parameter.h
 * @brief model::Parameter source file
 */
#ifndef SVZERODSOLVER_MODEL_PARAMETER_HPP_
#define SVZERODSOLVER_MODEL_PARAMETER_HPP_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "DOFHandler.h"

/**
 * @brief Model Parameter.
 *
 * This class handles constant parameters and time-dependent parameters that
 * need to be interpolated and periodically applied.
 *
 */
class Parameter {
 public:
  /**
   * @brief Construct a new Parameter object
   *
   * @param id Global ID of the parameter
   * @param value The value of the parameter
   */
  Parameter(int id, double value);

  /**
   * @brief Construct a new Parameter object
   *
   * @param id Global ID of the parameter
   * @param times Time steps corresponding to the time-dependent values
   * @param values Values corresponding to the time steps
   * @param periodic Is this parameter periodic with a cardiac cycle?
   */
  Parameter(int id, const std::vector<double>& times,
            const std::vector<double>& values, bool periodic = true);

  int id;                      ///< Global ID of the parameter
  std::vector<double> times;   ///< Time steps if parameter is time-dependent
  std::vector<double> values;  ///< Values if parameter is time-dependent
  double value;                ///< Value if parameter is constant
  double cycle_period;  ///< Cardiac cycle period corresponding to the time
                        ///< sequence
  int size;          ///< Size of the time series if parameter is time-dependent
  bool is_constant;  ///< Bool value indicating if the parameter is constant
  bool is_periodic;  ///< Bool value indicating if the parameter is periodic
                     ///< with the cardiac cycle

  /**
   * @brief Update the parameter
   *
   * @param value Value of the parameter
   */
  void update(double value);

  /**
   * @brief Update the parameter
   *
   * @param times Time steps corresponding to the values
   * @param values Values correspondong to the time steps
   */
  void update(const std::vector<double>& times,
              const std::vector<double>& values);

  /**
   * @brief Get the parameter value at the specified time.
   *
   * @param time Current time
   * @return Value at the time
   */
  double get(double time);

  /**
   * @brief Convert the parameter into a steady mean state.
   *
   */
  void to_steady();

  /**
   * @brief Convert the parameter back into an unsteady mean state.
   *
   */
  void to_unsteady();

 private:
  bool steady_converted = false;
};

/**
 * @brief Handles the properties of input parameters
 */
struct InputParameter {
  bool is_optional;    ///< Is this parameter optional?
  bool is_array;       ///< Is this parameter an array?
  bool is_number;      ///< Is this parameter a number?
  double default_val;  ///< Default value (if parameter is optional)

  /**
   * @brief Handles input parameters
   *
   * @param is_optional Is this parameter optional?
   * @param is_array Is this parameter an array?
   * @param is_number Is this parameter a number?
   * @param default_val Default value (if parameter is optional)
   */
  InputParameter(bool is_optional = false, bool is_array = false,
                 bool is_number = true, double default_val = 0.0)
      : is_optional(is_optional),
        is_array(is_array),
        is_number(is_number),
        default_val(default_val) {}
};

#endif  // SVZERODSOLVER_MODEL_PARAMETER_HPP_
