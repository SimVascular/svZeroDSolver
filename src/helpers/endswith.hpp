/**
 * @file endswith.hpp
 * @brief HELPERS::endswith source file
 */
#ifndef SVZERODSOLVER_HELPERS_ENDSWITH_HPP_
#define SVZERODSOLVER_HELPERS_ENDSWITH_HPP_

#include <string>

namespace HELPERS
{
    /**
     * @brief Check if a string ends with the letters of another string
     *
     * @param str The string to check for the specified suffix
     * @param suffix The suffix the string should be checked for
     * @return true if strings ends with the suffix, otherwise false
     */
    bool endswith(const std::string &str, const std::string &suffix)
    {
        if (str.length() >= suffix.length())
        {
            return (0 == str.compare(str.length() - suffix.length(), suffix.length(), suffix));
        }
        else
        {
            return false;
        }
    }
}

#endif // SVZERODSOLVER_HELPERS_ENDSWITH_HPP_