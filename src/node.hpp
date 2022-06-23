#ifndef SVZERODSOLVER_NODE_H_
#define SVZERODSOLVER_NODE_H_

#include "dofhandler.hpp"

class Block;

class Node
{
public:
    Node(Block &ele1, Block &ele2, std::string name);
    ~Node();
    std::string name;
    unsigned int flow_dof;
    unsigned int pres_dof;

    void setup_dofs(DOFHandler &dofhandler);
};

#include "block.hpp"

#endif // SVZERODSOLVER_NODE_H_