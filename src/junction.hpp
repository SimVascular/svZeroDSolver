#ifndef SVZERODSOLVER_JUNCTION_H_
#define SVZERODSOLVER_JUNCTION_H_

#include "block.hpp"

class Junction : public Block
{
protected:
    unsigned int num_internal_vars = 0;

public:
    struct Parameters : public Block::Parameters
    {
    };
    Junction(Parameters &params, std::string name);
    ~Junction();
    void setup_dofs(DOFHandler &dofhandler);
    void update_constant(System system);

private:
    Parameters *params;
};

#endif // SVZERODSOLVER_JUNCTION_H_