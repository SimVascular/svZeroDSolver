#ifndef SVZERODSOLVER_INTEGRATOR_H_
#define SVZERODSOLVER_INTEGRATOR_H_

#include <map>
#include <Eigen/Dense>
#include "system.hpp"
#include "model.hpp"

class State
{
public:
    Eigen::VectorXd y;
    Eigen::VectorXd ydot;
    State();
    ~State();
    State(const State &state);

    static State Zero(unsigned int n);
};

class Integrator
{
private:
    double alpha_m;
    double alpha_f;
    double gamma;
    double gamma_inv;
    double time_step_size;
    double time_step_size_inv;
    double y_dot_coeff;

    int size;
    Eigen::VectorXd y_af;
    Eigen::VectorXd ydot_am;
    Eigen::VectorXd dy;
    Eigen::MatrixXd lhs;
    Eigen::VectorXd rhs;

    System system;

public:
    Integrator(System &system, double time_step_size, double rho);
    ~Integrator();
    State step(State &state, double time, Model &model, unsigned int max_iter);
};

#endif // SVZERODSOLVER_INTEGRATOR_H_