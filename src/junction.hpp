#ifndef SVZERODSOLVER_JUNCTION_H_
#define SVZERODSOLVER_JUNCTION_H_

#include "block.hpp"

class Junction : public Block
{
public:
    struct Parameters : public Block::Parameters
    {
    };
    Junction(std::string name);
    ~Junction();
    void setup_dofs(DOFHandler &dofhandler);
    void update_constant(System &system);

    std::string name;

private:
    Parameters params;
    unsigned int num_inlets;
    unsigned int num_outlets;
};

#endif // SVZERODSOLVER_JUNCTION_H_