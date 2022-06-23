#include "dofhandler.hpp"
#include <iostream>

DOFHandler::DOFHandler()
{
    var_counter = 0;
    eqn_counter = 0;
}

DOFHandler::~DOFHandler()
{
}

unsigned int DOFHandler::size()
{
    return var_counter;
}

unsigned int DOFHandler::register_variable(std::string name)
{
    variables.push_back(name);
    return var_counter++;
}

unsigned int DOFHandler::register_equation()
{
    return eqn_counter++;
}
