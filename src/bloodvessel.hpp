#ifndef SVZERODSOLVER_BLOODVESSEL_H_
#define SVZERODSOLVER_BLOODVESSEL_H_

#include "block.hpp"

class BloodVessel : public Block
{
public:
    struct Parameters : public Block::Parameters
    {
        double R;                    // Poseuille resistance
        double C;                    // Capacitance
        double L;                    // Inductance
        double stenosis_coefficient; // Stenosis Coefficient
    };
    BloodVessel(double R, double C, double L, double stenosis_coefficient, std::string name);
    ~BloodVessel();
    void setup_dofs(DOFHandler &dofhandler);
    void update_constant(System &system);
    void update_solution(System &system, Eigen::VectorXd &y);

    std::string name;
    std::vector<Node *> inlet_nodes;
    std::vector<Node *> outlet_nodes;

private:
    Parameters params;
};

#endif // SVZERODSOLVER_BLOODVESSEL_H_