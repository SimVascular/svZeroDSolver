#include "bloodvessel.hpp"

BloodVessel::BloodVessel(BloodVessel::Parameters &params, std::string name) : Block(params, name)
{
    this->name = name;
    this->params = &params;
}

BloodVessel::~BloodVessel()
{
}

void BloodVessel::update_constant(System system)
{
    system.E(global_eqn_ids[0], global_var_ids[3]) = -params->L;
    system.E(global_eqn_ids[1], global_var_ids[4]) = -params->C;

    system.F(global_eqn_ids[0], global_var_ids[0]) = 1.0;
    system.F(global_eqn_ids[0], global_var_ids[1]) = -params->R;
    system.F(global_eqn_ids[0], global_var_ids[2]) = -1.0;
    system.F(global_eqn_ids[1], global_var_ids[1]) = 1.0;
    system.F(global_eqn_ids[1], global_var_ids[1]) = -1.0;
    system.F(global_eqn_ids[2], global_var_ids[0]) = 1.0;
    system.F(global_eqn_ids[2], global_var_ids[2]) = -params->R;
    system.F(global_eqn_ids[2], global_var_ids[4]) = -1.0;
}
void BloodVessel::update_solution(System system, Eigen::VectorXd &y)
{
    double q_in = abs(y[inlet_nodes[0]->flow_dof]);
    double fac1 = -params->stenosis_coefficient * q_in;
    double fac2 = fac1 - params->R;
    system.F(0, 1) = fac2;
    system.F(2, 1) = fac2;
    system.dF(0, 1) = fac1;
    system.dF(2, 1) = fac1;
}
