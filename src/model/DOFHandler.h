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
      output_variables;  ///< Names of the post-processing variables, in the
                         ///< order they occupy at the tail of the state

  /**
   * @brief Get the length of the state vector
   *
   * This counts every variable, the post-processing variables included.
   *
   * @return Length of the state vector
   */
  int size() const;

  /**
   * @brief Get the number of equations, which is the size of the non-linear
   * system that is actually solved
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
   * A post-processing variable is part of the state vector and is written to
   * the output, but it has no equation: it does not appear in any residual,
   * adds no row or column to the Jacobian and is never solved for. The block
   * that owns it computes it from the converged solution in
   * Block::update_outputs. Use this for anything that is reported but does
   * not feed back into the physics, such as the open/closed state of a
   * valve.
   *
   * Post-processing variables occupy the tail of the state vector, after
   * every solved variable, which is what lets the solved system stay a
   * leading square block. Model::finalize therefore registers them in a
   * second pass, once every block has registered its solved variables.
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
