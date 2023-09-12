
#include "helpers.h"

namespace helpers {

bool endswith(const std::string &str, const std::string &suffix) 
{
  if (str.length() >= suffix.length()) {
    return (0 == str.compare(str.length() - suffix.length(), suffix.length(),
                             suffix));
  } else {
    return false;
  }
}


bool startswith(const std::string &str, const std::string &prefix) 
{
  return str.size() >= prefix.size() &&
         str.compare(0, prefix.size(), prefix) == 0;
}

}

