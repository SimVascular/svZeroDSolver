/**
 * @file startswith.hpp
 * @brief HELPERS::startswith source file
 */
#ifndef SVZERODSOLVER_HELPERS_STARTSWITH_HPP_
#define SVZERODSOLVER_HELPERS_STARTSWITH_HPP_

#include <string>

namespace HELPERS
{
    /**
     * @brief Check if a string starts with the letters of another string
     *
     * @param str The string to check for the specified prefix
     * @param prefix The prefix the string should be checked for
     * @return true if strings starts with the prefix, otherwise false
     */
    bool startswith(const std::string &str, const std::string &prefix)
    {
        return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
    }
}

#endif // SVZERODSOLVER_HELPERS_STARTSWITH_HPP_