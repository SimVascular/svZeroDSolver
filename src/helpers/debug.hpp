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

#endif // SVZERODSOLVER_IO_CSVWRITER_HPP_