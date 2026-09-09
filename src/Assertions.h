// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file Assertions.h
 * @brief eigen_assert override for builds configured with ENABLE_ASSERTIONS
 */
#ifndef SVZERODSOLVER_ASSERTIONS_H_
#define SVZERODSOLVER_ASSERTIONS_H_

#include <stdexcept>
#include <string>

/**
 * @brief Report a violated Eigen assertion
 *
 * Eigen asserts via assert(), which aborts the process. Inside the pysvzerod
 * extension module that would tear down the Python interpreter, so throw
 * instead: pybind11 surfaces the exception and only the current call fails.
 *
 * @param condition Stringified condition that was violated
 * @param file Source file containing the assertion
 * @param line Line number of the assertion
 */
[[noreturn]] inline void svzerod_eigen_assertion_failed(const char* condition,
                                                        const char* file,
                                                        int line) {
  throw std::runtime_error("Eigen assertion failed: " + std::string(condition) +
                           " at " + std::string(file) + ":" +
                           std::to_string(line));
}

/**
 * @brief Eigen assertion that throws instead of aborting
 *
 * An expression, not a statement, because Eigen also uses eigen_assert inside
 * larger expressions.
 *
 * @param condition Condition that must hold
 */
#define eigen_assert(condition) \
  ((condition)                  \
       ? static_cast<void>(0)   \
       : svzerod_eigen_assertion_failed(#condition, __FILE__, __LINE__))

#endif  // SVZERODSOLVER_ASSERTIONS_H_
