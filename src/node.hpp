#ifndef SVZERODSOLVER_NODE_H_
#define SVZERODSOLVER_NODE_H_

#include "dofhandler.hpp"
#include <string>

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

Node::Node(std::string name)
{
    this->name = name;
}

Node::~Node()
{
}

void Node::setup_dofs(DOFHandler &dofhandler)
{
    flow_dof = dofhandler.register_variable("Q_" + name);
    pres_dof = dofhandler.register_variable("P_" + name);
}

#endif // SVZERODSOLVER_NODE_H_