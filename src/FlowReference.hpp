#ifndef SVZERODSOLVER_FLOWREFERENCE_H_
#define SVZERODSOLVER_FLOWREFERENCE_H_

#include "block.hpp"

class FlowReference : public Block
{
public:
    struct Parameters : public Block::Parameters
    {
        double Q; // Flow at timestep
    };
    FlowReference(Parameters &params, std::string name);
    ~FlowReference();
    void setup_dofs(DOFHandler &dofhandler);
    void update_constant(System system);
    void update_time(System system);

private:
    Parameters *params;
};

#endif // SVZERODSOLVER_FLOWREFERENCE_H_