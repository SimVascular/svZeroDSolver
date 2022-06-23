#ifndef SVZERODSOLVER_JUNCTION_H_
#define SVZERODSOLVER_JUNCTION_H_

#include "block.hpp"

class Junction : public Block
{
public:
    struct Parameters : public Block::Parameters
    {
    };
    Junction(Parameters &params, std::string name);
    ~Junction();
    void setup_dofs(DOFHandler &dofhandler);
    void update_constant(System system);

    std::string name;
    std::vector<Node *> inlet_nodes;
    std::vector<Node *> outlet_nodes;

private:
    Parameters *params;
};

#endif // SVZERODSOLVER_JUNCTION_H_