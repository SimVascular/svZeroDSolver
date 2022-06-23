#ifndef SVZERODSOLVER_RCRBLOCKWITHDISTALPRESSURE_H_
#define SVZERODSOLVER_RCRBLOCKWITHDISTALPRESSURE_H_

#include "block.hpp"

class RCRBlockWithDistalPressure : public Block
{
protected:
    unsigned int num_equations = 2;
    unsigned int num_internal_vars = 1;

public:
    struct Parameters
    {
        double Rp; // Proximal resistance
        double C;  // Capacitance
        double Rd; // Distal restistance
        double Pd; // Distal Pressure
    };
    void update_constant(System system);
    void update_time(System system);

private:
    Parameters *params;
};

#endif // SVZERODSOLVER_RCRBLOCKWITHDISTALPRESSURE_H_