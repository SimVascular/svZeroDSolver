// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "Block.h"

#include "Model.h"

std::string Block::get_name() const {
  return this->model->get_block_name(this->id);
}

void Block::update_vessel_type(VesselType type) { vessel_type = type; }

Block::~Block() {}

void Block::setup_params_(const std::vector<int>& param_ids) {
  this->global_param_ids = param_ids;
}

void Block::setup_dofs_(DOFHandler& dofhandler, int num_equations,
                        const std::list<std::string>& internal_var_names,
                        const std::list<std::string>& output_var_names) {
  // A block contributes one equation per variable it owns: each connection
  // contributes one of the two degrees of freedom of its node (the other
  // belongs to the block on the far side), plus one per internal variable.
  // Output variables are not solved for and have no equation, so they do not
  // count here. A block with too few connections would otherwise index past
  // the end of global_var_ids while assembling, which used to segfault.
  const int num_connections =
      static_cast<int>(inlet_nodes.size() + outlet_nodes.size());
  const int num_required =
      num_equations - static_cast<int>(internal_var_names.size());

  if (num_required < 0) {
    throw std::runtime_error(
        "Internal error: block '" + get_name() + "' declares " +
        std::to_string(num_equations) + " equations but " +
        std::to_string(internal_var_names.size()) + " internal variables.");
  }

  if (num_connections != num_required) {
    throw std::runtime_error(
        "Block '" + get_name() + "' needs " + std::to_string(num_required) +
        " connection" + (num_required == 1 ? "" : "s") +
        " (inlets plus outlets), but the model connects it " +
        std::to_string(num_connections) + " time" +
        (num_connections == 1 ? "" : "s") +
        ". Check that it is listed in a junction, or has a boundary "
        "condition, on every side.");
  }

  // Collect external DOFs from inlet nodes
  for (auto& inlet_node : inlet_nodes) {
    global_var_ids.push_back(inlet_node->pres_dof);
    global_var_ids.push_back(inlet_node->flow_dof);
  }

  // Collect external DOFs from outlet nodes
  for (auto& outlet_node : outlet_nodes) {
    global_var_ids.push_back(outlet_node->pres_dof);
    global_var_ids.push_back(outlet_node->flow_dof);
  }

  // Register internal variables of block
  for (auto& int_name : internal_var_names) {
    global_var_ids.push_back(
        dofhandler.register_variable(int_name + ":" + this->get_name()));
  }

  // Registered by Model::finalize in a second pass, so they land behind
  // every solved variable
  this->output_var_names = output_var_names;

  // Register equations of block
  for (int i = 0; i < num_equations; i++) {
    global_eqn_ids.push_back(dofhandler.register_equation(get_name()));
  }
}

void Block::setup_dofs(DOFHandler& dofhandler) {}

void Block::setup_model_dependent_params() {}

void Block::setup_initial_state_dependent_params(
    State initial_state, std::vector<double>& parameters) {}

void Block::update_constant(SparseSystem& system,
                            std::vector<double>& parameters) {}

void Block::update_time(SparseSystem& system, std::vector<double>& parameters) {
}

void Block::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {}

void Block::post_solve(Eigen::Matrix<double, Eigen::Dynamic, 1>& y) {}

void Block::update_outputs(
    Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& ydot) {}

void Block::accept_timestep(const Eigen::Matrix<double, Eigen::Dynamic, 1>& y) {
}

std::vector<double> Block::get_persistent_state() const { return {}; }

void Block::set_persistent_state(const std::vector<double>& state) {}

void Block::update_gradient(Eigen::SparseMatrix<double>& jacobian,
                            Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha,
                            std::vector<double>& y, std::vector<double>& dy) {
  throw std::runtime_error("Gradient calculation not implemented for block " +
                           get_name());
}

TripletsContributions Block::get_num_triplets() { return num_triplets; }
