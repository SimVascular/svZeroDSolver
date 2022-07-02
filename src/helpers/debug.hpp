/**
 * @file debug.hpp
 * @brief HELPERS::DEBUG source file
 */
#ifndef SVZERODSOLVER_HELPERS_DEBUG_HPP_
#define SVZERODSOLVER_HELPERS_DEBUG_HPP_

#include <iostream>

#ifdef DEBUG
#define DEBUG_MSG(str)                 \
    do                                 \
    {                                  \
        std::cout << str << std::endl; \
    } while (false)
#else
#define DEBUG_MSG(str) \
    do                 \
    {                  \
    } while (false)
#endif

#endif // SVZERODSOLVER_HELPERS_DEBUG_HPP_