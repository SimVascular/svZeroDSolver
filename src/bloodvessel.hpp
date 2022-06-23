#ifndef SVZERODSOLVER_BLOODVESSEL_H_
#define SVZERODSOLVER_BLOODVESSEL_H_

#include "block.hpp"

class BloodVessel : public Block
{
protected:
    unsigned int num_equations = 3;
    unsigned int num_internal_vars = 1;

public:
    struct Parameters
    {
        double R;                    // Poseuille resistance
        double C;                    // Capacitance
        double L;                    // Inductance
        double stenosis_coefficient; // Stenosis Coefficient
    };
    void update_constant(System system);
    void update_solution(System system, Eigen::VectorXd &y);

private:
    Parameters *params;
};

#endif // SVZERODSOLVER_BLOODVESSEL_H_