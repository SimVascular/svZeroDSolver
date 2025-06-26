// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "Node.h"

#include "Block.h"
#include "Model.h"

Node::Node(int id, const std::vector<Block *> &inlet_eles,
           const std::vector<Block *> &outlet_eles, Model *model) {
  this->id = id;
  this->inlet_eles = inlet_eles;
  this->outlet_eles = outlet_eles;
  this->model = model;

  for (auto &inlet_ele : inlet_eles) {
    inlet_ele->outlet_nodes.push_back(this);
  }

  for (auto &outlet_ele : outlet_eles) {
    outlet_ele->inlet_nodes.push_back(this);
  }
}

std::string Node::get_name() { return this->model->get_node_name(this->id); }

void Node::setup_dofs(DOFHandler &dofhandler) {
  flow_dof = dofhandler.register_variable("flow:" + get_name());
  pres_dof = dofhandler.register_variable("pressure:" + get_name());
}
