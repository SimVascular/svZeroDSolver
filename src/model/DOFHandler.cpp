
#include "DOFHandler.h"
#include <algorithm>  

namespace zd_model {

DOFHandler::DOFHandler() 
{
  var_counter = 0;
  eqn_counter = 0;
}

DOFHandler::~DOFHandler() {}

int DOFHandler::size() 
{ 
  return eqn_counter; 
}

int DOFHandler::get_num_equations() 
{ 
  return eqn_counter; 
}

int DOFHandler::get_num_variables() 
{ 
  return var_counter; 
}

int DOFHandler::register_variable(std::string name) 
{
  variables.push_back(name);
  variable_name_map.insert({name, var_counter});
  return var_counter++;
}

int DOFHandler::get_variable_index(std::string name) 
{
  return variable_name_map[name];
}

int DOFHandler::register_equation(std::string name) 
{
  equations.push_back(name);
  return eqn_counter++;
}

int DOFHandler::get_index(std::string_view& name) 
{
  auto it = std::find(variables.begin(), variables.end(), name);

  if (it != variables.end()) {
    return it - variables.begin();
  } else {
    throw std::runtime_error("No variable with that name");
  }
}

} 

