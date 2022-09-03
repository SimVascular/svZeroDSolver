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
/**
 * @file node.hpp
 * @brief MODEL::Node source file
 */
#ifndef SVZERODSOLVER_MODEL_NODE_HPP_
#define SVZERODSOLVER_MODEL_NODE_HPP_

#include <string>

#include "dofhandler.hpp"

namespace MODEL {

/**
 * @brief Node
 *
 * Nodes connect two blocks with each other. Each node corresponds to a
 * flow and pressure value of the system.
 */
class Node {
 public:
  /**
   * @brief Construct a new Node object
   *
   * @param name Name
   */
  Node(std::string name);

  /**
   * @brief Destroy the Node object
   *
   */
  ~Node();

  std::string name;       ///< Name
  unsigned int flow_dof;  ///< Global flow degree-of-freedom of the node
  unsigned int pres_dof;  ///< Global pressure degree-of-freedom of the node

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set \ref global_var_ids and \ref global_eqn_ids of the element based on the
   * number of equations and the number of internal variables of the
   * element.
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler &dofhandler);
};

Node::Node(std::string name) { this->name = name; }

Node::~Node() {}

void Node::setup_dofs(DOFHandler &dofhandler) {
  flow_dof = dofhandler.register_variable("flow:" + name);
  pres_dof = dofhandler.register_variable("pressure:" + name);
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_NODE_HPP_