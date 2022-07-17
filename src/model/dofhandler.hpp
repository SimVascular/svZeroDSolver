/**
 * @file dofhandler.hpp
 * @brief MODEL::DOFHandler source file
 */
#ifndef SVZERODSOLVER_MODEL_DOFHANDLER_HPP_
#define SVZERODSOLVER_MODEL_DOFHANDLER_HPP_

#include <string>
#include <vector>

namespace MODEL {

/**
 * @brief Degree-of-freedom handler.
 *
 * This class handles degrees-of-freedom for model variables and
 * equations. It assigns each element with row and column indices which it
 * can use to assemble it's local contributions into the global system.
 */
class DOFHandler {
 private:
  unsigned int var_counter;  ///< Variable counter
  unsigned int eqn_counter;  ///< Equation counter

 public:
  std::vector<std::string>
      variables;  ///< Variable names corresponding to the variable indices

  /**
   * @brief Construct a new DOFHandler object
   *
   */
  DOFHandler();

  /**
   * @brief Destroy the DOFHandler object
   *
   */
  ~DOFHandler();

  /**
   * @brief Get the size of the system
   *
   * @return Size of the system
   */
  unsigned int size();

  /**
   * @brief Register a new variable at the DOFHandler.
   *
   * @param name Name of the variable
   * @return Global index of the variable
   */
  unsigned int register_variable(std::string name);

  /**
   * @brief Register a new equation at the DOFHandler
   *
   * @return Global index of the equation
   */
  unsigned int register_equation();
};

DOFHandler::DOFHandler() {
  var_counter = 0;
  eqn_counter = 0;
}

DOFHandler::~DOFHandler() {}

unsigned int DOFHandler::size() { return var_counter; }

unsigned int DOFHandler::register_variable(std::string name) {
  variables.push_back(name);
  return var_counter++;
}

unsigned int DOFHandler::register_equation() { return eqn_counter++; }

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_DOFHANDLER_HPP_