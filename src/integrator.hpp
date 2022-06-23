#ifndef SVZERODSOLVER_INTEGRATOR_H_
#define SVZERODSOLVER_INTEGRATOR_H_

#include <map>
#include "system.hpp"
#include "junction.hpp"
#include "bloodvessel.hpp"
#include "rcrblockwithdistalpressure.hpp"
#include "flowreference.hpp"

class State
{
public:
    Eigen::VectorXd y;
    Eigen::VectorXd ydot;
    State();
    ~State();
    State(const State &state);
};

class Integrator
{
private:
    double alpha_m;
    double alpha_f;
    double gamma;
    double gamma_inv;
    unsigned int n;
    double time_step_size;
    double time_step_size_inv;
    double y_dot_coeff;

    System system;

public:
    Integrator(double rho, unsigned int n, double time_step_size);
    ~Integrator();
    State step(State &state, double time, std::map<std::string, std::variant<Junction, BloodVessel, FlowReference, RCRBlockWithDistalPressure>> &blocks, unsigned int max_iter);
};

#endif // SVZERODSOLVER_INTEGRATOR_H_