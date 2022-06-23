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

#endif // SVZERODSOLVER_DOFHANDLER_H_