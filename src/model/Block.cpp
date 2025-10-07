// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "Block.h"

#include "Model.h"

std::string Block::get_name() { return this->model->get_block_name(this->id); }

void Block::update_vessel_type(VesselType type) { vessel_type = type; }

Block::~Block() {}

void Block::setup_params_(const std::vector<int> &param_ids) {
  this->global_param_ids = param_ids;
}

void Block::setup_dofs_(DOFHandler &dofhandler, int num_equations,
                        const std::list<std::string> &internal_var_names) {
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
  for (int i = 0; i < num_equations; i++) {
    global_eqn_ids.push_back(dofhandler.register_equation(get_name()));
  }
}

void Block::setup_dofs(DOFHandler &dofhandler) {}

void Block::setup_model_dependent_params() {}

void Block::setup_initial_state_dependent_params(
    State initial_state, std::vector<double> &parameters) {}

void Block::update_constant(SparseSystem &system,
                            std::vector<double> &parameters) {}

void Block::update_time(SparseSystem &system, std::vector<double> &parameters) {
}

void Block::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {}

void Block::post_solve(Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {}

void Block::update_gradient(Eigen::SparseMatrix<double> &jacobian,
                            Eigen::Matrix<double, Eigen::Dynamic, 1> &residual,
                            Eigen::Matrix<double, Eigen::Dynamic, 1> &alpha,
                            std::vector<double> &y, std::vector<double> &dy) {
  throw std::runtime_error("Gradient calculation not implemented for block " +
                           get_name());
}

TripletsContributions Block::get_num_triplets() { return num_triplets; }
