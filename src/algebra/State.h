// This file is part of svZeroDSolver licensed under the MIT License
// See the LICENSE.md file for license information
/**
 * @file State.h
 * @brief State source file
 */
#ifndef SVZERODSOLVER_ALGEBRA_STATE_HPP_
#define SVZERODSOLVER_ALGEBRA_STATE_HPP_

#include <Eigen/Core>

/**
 * @brief State of the system.
 *
 * Stores the current state of a system, i.e. the current value and
 * derivate of all variables.
 *
 */
class State {
 public:
  Eigen::Matrix<double, Eigen::Dynamic, 1>
      y;  ///< Vector of solution quantities
  Eigen::Matrix<double, Eigen::Dynamic, 1> ydot;  ///< Derivate of \ref y

  /**
   * @brief Construct a new State object
   *
   */
  State();

  /**
   * @brief Construct a new State object
   *
   * @param n Size of the state
   */
  State(int n);

  /**
   * @brief Destroy the State object
   *
   */
  ~State();

  /**
   * @brief Copy a State object
   *
   * @param state
   */
  State(const State &state);

  /**
   * @brief Construct a new State object and initilaize with all zeros.
   *
   * @param n Size of the state
   * @return New state initialized with all zeros
   */
  static State Zero(int n);
};

#endif  // SVZERODSOLVER_ALGEBRA_STATE_HPP_
