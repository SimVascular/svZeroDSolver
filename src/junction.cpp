#include "junction.hpp"

void Junction::setup_dofs(DOFHandler &dofhandler)
{
    // Derive number of inlets and outlets
    unsigned int num_inlets = this->inlet_nodes.size();
    unsigned int num_outlets = this->outlet_nodes.size();

    // Set number of equations of a junction block based on number of
    // inlets/outlets. Must be set before calling parent constructor
    this->num_equations = num_inlets + num_outlets;
    Block::setup_dofs(dofhandler);
}

void Junction::update_constant(System system)
{
    for (size_t i = 0; i < (num_equations - 1); i++)
    {
        system.F(global_eqn_ids[i], global_var_ids[0]) = 1.0;
        system.F(global_eqn_ids[i], global_var_ids[2 * i + 2]) = -1.0;
    }
    unsigned int num_inlets = this->inlet_nodes.size();
    unsigned int num_outlets = this->outlet_nodes.size();
    for (size_t i = 0; i < num_inlets; i++)
    {
        system.F(global_eqn_ids[num_equations - 1], global_var_ids[i]) = 1.0;
    }
    for (size_t i = num_inlets; i < num_inlets + num_outlets; i++)
    {
        system.F(global_eqn_ids[num_equations - 1], global_var_ids[i]) = -1.0;
    }
}