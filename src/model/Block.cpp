
#include "Block.h"
#include "Model.h"

namespace zd_model {

Block::Block(int id, const std::vector<int> &param_ids, Model *model) 
{
  this->id = id;
  this->global_param_ids = param_ids;
  this->model = model;
}

std::string Block::get_name() 
{
  return this->model->get_block_name(this->id);
}

Block::~Block() {}

void Block::setup_dofs_(DOFHandler &dofhandler, unsigned int num_equations,
                           std::list<std::string> internal_var_names) 
{
  // Collect external DOFs from inlet nodes
  for (auto &inlet_node : inlet_nodes) {
    global_var_ids.push_back(inlet_node->pres_dof);
    global_var_ids.push_back(inlet_node->flow_dof);
  }

  // Collect external DOFs from outlet nodes
  for (auto &outlet_node : outlet_nodes) {
    global_var_ids.push_back(outlet_node->pres_dof);
    global_var_ids.push_back(outlet_node->flow_dof);
  }

  // Register internal variables of block
  for (auto &int_name : internal_var_names) {
    global_var_ids.push_back(
        dofhandler.register_variable(int_name + ":" + this->get_name()));
  }

  // Register equations of block
  for (unsigned int i = 0; i < num_equations; i++) {
    global_eqn_ids.push_back(dofhandler.register_equation(get_name()));
  }
}

void Block::setup_dofs(DOFHandler &dofhandler) {}

void Block::setup_model_dependent_params() {}

void Block::update_constant(algebra::SparseSystem &system, std::vector<double> &parameters) {}

void Block::update_time(algebra::SparseSystem& system,
                           std::vector<double> &parameters) {}

void Block::update_solution(algebra::SparseSystem& system, std::vector<double> &parameters,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &y, Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {}

void Block::update_gradient(Eigen::SparseMatrix<double> &jacobian, 
    Eigen::Matrix<double, Eigen::Dynamic, 1> &residual, Eigen::Matrix<double, Eigen::Dynamic, 1> &alpha,
    std::vector<double> &y, std::vector<double> &dy) 
{
  throw std::runtime_error("Gradient calculation not implemented for block " +
                           get_name());
}

std::map<std::string, int> Block::get_num_triplets() {
  return num_triplets;
}

};  

