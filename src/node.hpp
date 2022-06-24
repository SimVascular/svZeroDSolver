#ifndef SVZERODSOLVER_NODE_H_
#define SVZERODSOLVER_NODE_H_

#include "dofhandler.hpp"

class Node
{
public:
    Node(std::string name);
    ~Node();
    std::string name;
    unsigned int flow_dof;
    unsigned int pres_dof;

    void setup_dofs(DOFHandler &dofhandler);
};

#endif // SVZERODSOLVER_NODE_H_