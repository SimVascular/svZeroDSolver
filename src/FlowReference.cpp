#include "flowreference.hpp"

#include <iostream>

FlowReference::FlowReference(TimeDependentParameter Q, std::string name) : Block(name)
{
    this->name = name;
    this->params.Q = Q;
}

FlowReference::~FlowReference()
{
}

void FlowReference::setup_dofs(DOFHandler &dofhandler)
{
    Block::setup_dofs_(dofhandler, 1, 0);
}

void FlowReference::update_constant(System &system)
{
    std::cout << "FlowReference" << std::endl;
    system.F(global_eqn_ids[0], global_var_ids[1]) = 1.0;
}
void FlowReference::update_time(System &system, double time)
{
    system.C(global_eqn_ids[0]) = -params.Q.get(time);
}