// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "DOFHandler.h"

#include <algorithm>  
#include <stdexcept>

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

