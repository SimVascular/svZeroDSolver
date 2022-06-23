#include "block.hpp"

Block::Block(Block::Parameters &params, std::string name)
{
    this->name = name;
    this->params = &params;
}

Block::~Block()
{
}

void Block::setup_dofs(DOFHandler &dofhandler)
{
    // Collect external DOFs from inlet and outlet nodes
    for (auto it = inlet_nodes.begin(); it != inlet_nodes.end(); it++)
    {
        global_var_ids.push_back((*it)->pres_dof);
        global_var_ids.push_back((*it)->flow_dof);
    }
    for (auto it = outlet_nodes.begin(); it != outlet_nodes.end(); it++)
    {
        global_var_ids.push_back((*it)->pres_dof);
        global_var_ids.push_back((*it)->flow_dof);
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

// Bloodvessel

// FlowRef BC

// RCR BC
