// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Block.h"

#include "Model.h"

Block::Block(int id, const std::vector<int> &param_ids, Model *model) {
  this->id = id;
  this->global_param_ids = param_ids;
  this->model = model;
}

std::string Block::get_name() { return this->model->get_block_name(this->id); }

Block::~Block() {}

void Block::setup_dofs_(DOFHandler &dofhandler, unsigned int num_equations,
                        std::list<std::string> internal_var_names) {
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

void Block::update_constant(SparseSystem &system,
                            std::vector<double> &parameters) {}

void Block::update_time(SparseSystem &system,
                        std::vector<double> &parameters) {}

void Block::update_solution(SparseSystem &system,
                            std::vector<double> &parameters,
                            Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                            Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {}

void Block::update_gradient(Eigen::SparseMatrix<double> &jacobian,
                            Eigen::Matrix<double, Eigen::Dynamic, 1> &residual,
                            Eigen::Matrix<double, Eigen::Dynamic, 1> &alpha,
                            std::vector<double> &y, std::vector<double> &dy) {
  throw std::runtime_error("Gradient calculation not implemented for block " +
                           get_name());
}

std::map<std::string, int> Block::get_num_triplets() { return num_triplets; }

