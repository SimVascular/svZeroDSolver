#include "rcrblockwithdistalpressure.hpp"

RCRBlockWithDistalPressure::RCRBlockWithDistalPressure(double Rp, double C, double Rd, double Pd, std::string name) : Block(name)
{
    this->name = name;
    this->params.Rp = Rp;
    this->params.C = C;
    this->params.Rd = Rd;
    this->params.Pd = Pd;
}

RCRBlockWithDistalPressure::~RCRBlockWithDistalPressure()
{
}

void RCRBlockWithDistalPressure::setup_dofs(DOFHandler &dofhandler)
{
    Block::setup_dofs_(dofhandler, 2, 1);
}

void RCRBlockWithDistalPressure::update_constant(System &system)
{
    system.F(global_eqn_ids[0], global_var_ids[0]) = 1.0;
    system.F(global_eqn_ids[0], global_var_ids[2]) = -1.0;
    system.F(global_eqn_ids[1], global_var_ids[2]) = -1.0;
}
void RCRBlockWithDistalPressure::update_time(System &system, double time)
{
    system.E(global_eqn_ids[1], global_var_ids[2]) = -params.Rd * params.C;
    system.F(global_eqn_ids[0], global_var_ids[1]) = -params.Rp;
    system.F(global_eqn_ids[1], global_var_ids[1]) = params.Rd;
    system.C(global_eqn_ids[1]) = params.Pd;
}

void RCRBlockWithDistalPressure::to_steady()
{
    params.C = 0.0;
}
