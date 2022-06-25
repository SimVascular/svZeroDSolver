#ifndef SVZERODSOLVER_DOFHANDLER_H_
#define SVZERODSOLVER_DOFHANDLER_H_

#include <vector>
#include <string>

class DOFHandler
{
private:
    unsigned int var_counter;
    unsigned int eqn_counter;

public:
    std::vector<std::string> variables;
    DOFHandler();
    ~DOFHandler();

    unsigned int size();

    unsigned int register_variable(std::string name);

    unsigned int register_equation();
};

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

#endif // SVZERODSOLVER_DOFHANDLER_H_