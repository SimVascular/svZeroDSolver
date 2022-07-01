#ifndef SVZERODSOLVER_HELPERS_STARTSWITH_HPP_
#define SVZERODSOLVER_HELPERS_STARTSWITH_HPP_

#include <string>

namespace HELPERS
{
    bool startswith(const std::string &str, const std::string &prefix)
    {
        return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
    }
}

#endif // SVZERODSOLVER_HELPERS_STARTSWITH_HPP_