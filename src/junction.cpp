#include "junction.hpp"

#include <iostream>

Junction::Junction(std::string name) : Block(name)
{
    this->name = name;
}

Junction::~Junction()
{
}

void Junction::setup_dofs(DOFHandler &dofhandler)
{
    // Set number of equations of a junction block based on number of
    // inlets/outlets. Must be set before calling parent constructor
    num_inlets = inlet_nodes.size();
    num_outlets = outlet_nodes.size();
    Block::setup_dofs_(dofhandler, num_inlets + num_outlets, 0);
}

void Junction::update_constant(System &system)
{
    std::cout << "Junction" << std::endl;
    for (size_t i = 0; i < (num_inlets + num_outlets - 1); i++)
    {
        std::cout << 1 << std::endl;
        std::cout << "Num rows " << system.F.rows() << std::endl;
        std::cout << "Row " << global_eqn_ids[i] << std::endl;
        std::cout << "Col " << global_var_ids[0] << std::endl;
        system.F(global_eqn_ids[i], global_var_ids[0]) = 1.0;
        std::cout << 2 << std::endl;
        system.F(global_eqn_ids[i], global_var_ids[2 * i + 2]) = -1.0;
    }
    for (size_t i = 0; i < num_inlets; i++)
    {
        std::cout << 3 << std::endl;
        system.F(global_eqn_ids[num_inlets + num_outlets - 1], global_var_ids[i]) = 1.0;
    }
    for (size_t i = num_inlets; i < num_inlets + num_outlets; i++)
    {
        std::cout << 4 << std::endl;
        system.F(global_eqn_ids[num_inlets + num_outlets - 1], global_var_ids[i]) = -1.0;
    }
}