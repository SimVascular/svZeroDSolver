
#include <string>

#include "node.hpp"

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