
#include <string>

#include "node.hpp"

Node::Node(Block &ele1, Block &ele2, std::string name)
{
    ele1.outlet_nodes.push_back(this);
    ele2.inlet_nodes.push_back(this);
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