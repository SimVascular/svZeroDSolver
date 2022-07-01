#ifndef SVZERODSOLVER_MODEL_NODE_HPP_
#define SVZERODSOLVER_MODEL_NODE_HPP_

#include <string>

#include "dofhandler.hpp"

namespace MODEL
{

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

} // namespace MODEL

#endif // SVZERODSOLVER_MODEL_NODE_HPP_