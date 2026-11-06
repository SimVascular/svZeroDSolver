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
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "DOFHandler.h"

/**
 * @brief Model Parameter.
 *
 * This class handles three kinds of parameters:
 *
 * * constant parameters,
 * * time-dependent parameters given as a time series that is linearly
 *   interpolated and (optionally) repeated every cardiac cycle, and
 * * time-dependent parameters given as a mathematical expression of the
 *   simulation time `t` (e.g. `"2*pi*cos(2*pi*t)"`), compiled once with
 *   ExprTk and evaluated at every time step.
 *
 * An expression is evaluated at the absolute simulation time, not the time
 * within the current cardiac cycle, so it must itself be periodic with the
 * cardiac cycle if the simulation is meant to converge to a periodic state.
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

  /**
   * @brief Construct a new Parameter object
   *
   * @param id Global ID of the parameter
   * @param expression_string Expression of the simulation time `t`
   */
  Parameter(int id, const std::string& expression_string);

  /**
   * @brief Destructor
   */
  ~Parameter();

  /**
   * @brief Move constructor
   */
  Parameter(Parameter&&) noexcept;

  /**
   * @brief Move assignment operator
   * @return Reference to this object
   */
  Parameter& operator=(Parameter&&) noexcept;

  int id;                      ///< Global ID of the parameter
  std::vector<double> times;   ///< Time steps if parameter is time-dependent
  std::vector<double> values;  ///< Values if parameter is time-dependent
  double value = 0.0;          ///< Value if parameter is constant
  double cycle_period = 0.0;   ///< Cardiac cycle period corresponding to the
                               ///< time sequence
  int size = 0;  ///< Size of the time series if parameter is time-dependent
  bool is_constant = true;   ///< Bool value indicating if the parameter is
                             ///< constant
  bool is_periodic = true;   ///< Bool value indicating if the parameter is
                             ///< periodic with the cardiac cycle
  bool is_function = false;  ///< Is the parameter an expression of time?

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
   * @brief Update the parameter with an expression of the simulation time `t`
   *
   * Throws std::runtime_error with the ExprTk diagnostics if the expression
   * does not compile.
   *
   * @param update_string Expression of the simulation time `t`
   */
  void update(const std::string& update_string);

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
   * A time series is replaced by the mean of its values. An expression is
   * replaced by its time average over one cardiac cycle, because its own
   * period is unknown.
   *
   * @param cardiac_cycle_period Period to average expressions over
   */
  void to_steady(double cardiac_cycle_period);

  /**
   * @brief Convert the parameter back into an unsteady mean state.
   *
   */
  void to_unsteady();

 private:
  bool steady_converted = false;
  struct ExprtkImpl;                         ///< Forward declared PIMPL struct
  std::unique_ptr<ExprtkImpl> exprtk_impl_;  ///< PIMPL for exprtk types
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
