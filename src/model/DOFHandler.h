// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
