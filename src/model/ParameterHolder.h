// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file ParameterHolder.h
 * @brief Common base for chamber sub-components configured from input
 * parameters
 */

#ifndef SVZERODSOLVER_MODEL_PARAMETERHOLDER_HPP_
#define SVZERODSOLVER_MODEL_PARAMETERHOLDER_HPP_

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Parameter.h"

/**
 * @brief Base class for chamber sub-components (activation functions,
 * materials, active stress models) that are configured from scalar/vector
 * input parameters read from the JSON input file.
 */
class ParameterHolder {
 public:
  /**
   * @brief Properties of the input parameters for this component
   * [(name, InputParameter), ...]
   */
  const std::vector<std::pair<std::string, InputParameter>>
      input_param_properties;

  /**
   * @brief Construct a ParameterHolder
   *
   * @param props Parameter name/spec pairs for this component, used to
   * initialize params_/params_vec_ with default values
   */
  explicit ParameterHolder(
      const std::vector<std::pair<std::string, InputParameter>>& props);

  virtual ~ParameterHolder() = default;

  /**
   * @brief Set a scalar parameter value by name.
   *
   * Calling function must validate the parameter name and value.
   *
   * @param name Parameter name
   * @param value Parameter value
   */
  void set_param(const std::string& name, double value);

  /**
   * @brief Set a vector-valued parameter value by name.
   *
   * Calling function must validate the parameter name and value.
   *
   * @param name Parameter name
   * @param value Parameter value
   */
  void set_param_vector(const std::string& name,
                        const std::vector<double>& value);

  /**
   * @brief Called once after all parameters have been read from the input
   * file (e.g. by the loader).
   *
   * Default no-op; subclasses override to validate the parameters they were
   * given or precompute derived values (e.g. normalization constants).
   */
  virtual void finalize() {}

 protected:
  /**
   * @brief Map of scalar parameter names to their values
   */
  std::map<std::string, double> params_;

  /**
   * @brief Map of vector-valued parameter names to their values
   */
  std::map<std::string, std::vector<double>> params_vec_;
};

#endif  // SVZERODSOLVER_MODEL_PARAMETERHOLDER_HPP_
