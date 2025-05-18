// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "DOFHandler.h"

#include <algorithm>
#include <stdexcept>

int DOFHandler::size() const { return eqn_counter; }

int DOFHandler::get_num_equations() const { return eqn_counter; }

int DOFHandler::get_num_variables() const { return var_counter; }

int DOFHandler::register_variable(const std::string& name) {
  variables.push_back(name);
  variable_name_map.insert({name, var_counter});
  return var_counter++;
}

int DOFHandler::get_variable_index(const std::string& name) const {
  try {
    return variable_name_map.at(name);
  } catch (...) {
    std::string error_msg = "ERROR: Variable name '" + name + "' not found.";
    throw std::runtime_error(error_msg);
  }
}

int DOFHandler::register_equation(const std::string& name) {
  equations.push_back(name);
  return eqn_counter++;
}

int DOFHandler::get_index(const std::string_view& name) const {
  auto it = std::find(variables.begin(), variables.end(), name);

  if (it != variables.end()) {
    return it - variables.begin();
  } else {
    throw std::runtime_error("No variable with that name");
  }
}
