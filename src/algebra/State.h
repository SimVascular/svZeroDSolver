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
