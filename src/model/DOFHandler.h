// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file DofHandler.h
 * @brief model::DOFHandler source file
 */
#ifndef SVZERODSOLVER_MODEL_DOFHANDLER_HPP_
#define SVZERODSOLVER_MODEL_DOFHANDLER_HPP_

#include <map>
#include <string>
#include <vector>

/**
 * @brief Degree-of-freedom handler.
 *
 * This class handles degrees-of-freedom for model variables and
 * equations. It assigns each element with row and column indices which it
 * can use to assemble it's local contributions into the global system.
 */
class DOFHandler {
 private:
  int var_counter{0};  ///< Variable counter
  int eqn_counter{0};  ///< Equation counter

 public:
  std::vector<std::string>
      variables;  ///< Variable names corresponding to the variable indices
  std::map<std::string, int>
      variable_name_map;  ///< Map between variable name and index
  std::vector<std::string>
      equations;  ///< Equation names corresponding to the equation indices
  std::vector<bool> equation_is_output;  ///< Marks equations that only define
                                         ///< a diagnostic output and do not
                                         ///< feed back into the physics

  /**
   * @brief Get the size of the system
   *
   * @return Size of the system
   */
  int size() const;

  /**
   * @brief Get the number of equations
   *
   * @return int Number of equations
   */
  int get_num_equations() const;

  /**
   * @brief Get the number of variables
   *
   * @return int Number of variables
   */
  int get_num_variables() const;

  /**
   * @brief Register a new variable at the DOFHandler.
   *
   * @param name Name of the variable
   * @return Global index of the variable
   */
  int register_variable(const std::string& name);

  /**
   * @brief Get the index of a variable by its name
   *
   * @param name Name of the variable
   * @return int Name of the variable
   */
  int get_variable_index(const std::string& name) const;

  /**
   * @brief Register a new equation at the DOFHandler
   *
   * @param name Name of the equation
   * @param is_output True if the equation only defines a diagnostic output.
   * Such equations are assembled and solved like any other, but their
   * residual does not take part in the convergence check of the non-linear
   * iteration. This matters for outputs that are discontinuous functions of
   * the solution, such as the open/closed state of an ideal valve: the
   * indicator flips by one whenever the pressure difference across the valve
   * changes sign, which would otherwise keep a converged time step from
   * terminating.
   * @return Global index of the equation
   */
  int register_equation(const std::string& name, bool is_output = false);

  /**
   * @brief Get the index of a variable
   *
   * @param name Name of the variable
   * @return Index of variable with given name
   */
  int get_index(const std::string_view& name) const;
};

#endif  // SVZERODSOLVER_MODEL_DOFHANDLER_HPP_
