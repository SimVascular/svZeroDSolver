// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file Debug.h
 * @brief DEBUG_MSG source file
 */
#ifndef SVZERODSOLVER_HELPERS_DEBUG_HPP_
#define SVZERODSOLVER_HELPERS_DEBUG_HPP_

#include <iostream>

/**
 * @brief DEBUG_MSG Macro to print debug messages for debug build
 */
#ifndef NDEBUG
#define DEBUG_MSG(str)                                   \
  do {                                                   \
    std::cout << "[DEBUG MESSAGE] " << str << std::endl; \
  } while (false)
#else
#define DEBUG_MSG(str) \
  do {                 \
  } while (false)
#endif

#endif  // SVZERODSOLVER_HELPERS_DEBUG_HPP_
