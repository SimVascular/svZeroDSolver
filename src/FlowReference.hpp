#ifndef SVZERODSOLVER_FLOWREFERENCE_H_
#define SVZERODSOLVER_FLOWREFERENCE_H_

#include "block.hpp"

class FlowReference : public Block
{
protected:
    unsigned int num_equations = 1;

public:
    struct Parameters
    {
        double Q; // Flow at timestep
    };
    void update_constant(System system);
    void update_time(System system);

private:
    Parameters *params;
};

#endif // SVZERODSOLVER_FLOWREFERENCE_H_