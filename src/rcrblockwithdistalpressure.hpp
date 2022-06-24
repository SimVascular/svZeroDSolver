#ifndef SVZERODSOLVER_RCRBLOCKWITHDISTALPRESSURE_H_
#define SVZERODSOLVER_RCRBLOCKWITHDISTALPRESSURE_H_

#include "block.hpp"

class RCRBlockWithDistalPressure : public Block
{

public:
    struct Parameters : public Block::Parameters
    {
        double Rp; // Proximal resistance
        double C;  // Capacitance
        double Rd; // Distal restistance
        double Pd; // Distal Pressure
    };
    RCRBlockWithDistalPressure(double Rp, double C, double Rd, double Pd, std::string name);
    ~RCRBlockWithDistalPressure();
    void setup_dofs(DOFHandler &dofhandler);
    void update_constant(System &system);
    void update_time(System &system, double time);

private:
    Parameters params;
};

#endif // SVZERODSOLVER_RCRBLOCKWITHDISTALPRESSURE_H_