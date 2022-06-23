#include "rcrblockwithdistalpressure.hpp"

RCRBlockWithDistalPressure::RCRBlockWithDistalPressure(RCRBlockWithDistalPressure::Parameters &params, std::string name) : Block(params, name)
{
    this->name = name;
    this->params = &params;
}

RCRBlockWithDistalPressure::~RCRBlockWithDistalPressure()
{
}

void RCRBlockWithDistalPressure::setup_dofs(DOFHandler &dofhandler)
{
    Block::setup_dofs_(dofhandler, 2, 1);
}

void RCRBlockWithDistalPressure::update_constant(System system)
{
    system.F(global_eqn_ids[0], global_var_ids[0]) = 1.0;
    system.F(global_eqn_ids[0], global_var_ids[2]) = -1.0;
    system.F(global_eqn_ids[1], global_var_ids[2]) = -1.0;
}
void RCRBlockWithDistalPressure::update_time(System system)
{
    system.E(global_eqn_ids[1], global_var_ids[2]) = -params->Rd * params->C;
    system.F(global_eqn_ids[0], global_var_ids[1]) = -params->Rp;
    system.F(global_eqn_ids[1], global_var_ids[1]) = params->Rd;
    system.C(global_eqn_ids[1]) = params->Pd;
}