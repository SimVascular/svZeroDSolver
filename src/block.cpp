#include "block.hpp"
#include <iostream>

Block::Block()
{
}

Block::Block(Block::Parameters &params, std::string name)
{
    this->name = name;
    this->params = &params;
}

Block::~Block()
{
}

void Block::setup_dofs_(DOFHandler &dofhandler, unsigned int num_equations, unsigned int num_internal_vars)
{
    // Collect external DOFs from inlet and outlet nodes
    for (size_t i = 0; i < inlet_nodes.size(); i++)
    {
        global_var_ids.push_back(inlet_nodes[i]->pres_dof);
        global_var_ids.push_back(inlet_nodes[i]->flow_dof);
    }
    for (size_t i = 0; i < outlet_nodes.size(); i++)
    {
        global_var_ids.push_back(outlet_nodes[i]->pres_dof);
        global_var_ids.push_back(outlet_nodes[i]->flow_dof);
    }

    // Register internal variables of block
    for (unsigned int i = 0; i < num_internal_vars; i++)
    {
        global_var_ids.push_back(dofhandler.register_variable("var_" + std::to_string(i) + "_" + name));
    }

    // Register equations of block
    for (unsigned int i = 0; i < num_equations; i++)
    {
        global_eqn_ids.push_back(dofhandler.register_equation());
    }
}
void Block::update_constant(System system)
{
}
void Block::update_time(System system, double time)
{
}
void Block::update_solution(System system, Eigen::VectorXd &y)
{
}
