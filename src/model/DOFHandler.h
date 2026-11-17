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
  std::vector<std::string>
      output_variables;  ///< Post-processing variables, at the state tail

  /**
   * @brief Get the number of equations, the size of the solved system
   *
   * @return int Number of equations
   */
  int get_num_equations() const;

  /**
   * @brief Get the number of post-processing variables
   *
   * @return int Number of post-processing variables
   */
  int get_num_output_variables() const;

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
   * @brief Register a new post-processing variable at the DOFHandler.
   *
   * Such a variable is in the state vector and in the output, but has no
   * equation: no row, no column, never solved for. The owning block computes
   * it in Block::update_outputs. They occupy the tail of the state vector so
   * the solved system stays a leading square block, which Model::finalize
   * arranges with a second registration pass.
   *
   * @param name Name of the variable
   * @return Global index of the variable
   */
  int register_output_variable(const std::string& name);

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
   * @return Global index of the equation
   */
  int register_equation(const std::string& name);

  /**
   * @brief Get the index of a variable
   *
   * @param name Name of the variable
   * @return Index of variable with given name
   */
  int get_index(const std::string_view& name) const;
};

#endif  // SVZERODSOLVER_MODEL_DOFHANDLER_HPP_
