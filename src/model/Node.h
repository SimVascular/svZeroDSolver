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
 * @file Node.h
 * @brief model::Node source file
 */
#ifndef SVZERODSOLVER_MODEL_NODE_HPP_
#define SVZERODSOLVER_MODEL_NODE_HPP_

#include <string>
#include <vector>

#include "BlockType.h"
#include "DOFHandler.h"

class Block;
class Model;

/**
 * @brief Node
 *
 * Nodes connect two blocks with each other. Each node corresponds to a
 * flow and pressure value of the system.
 *
 */
class Node {
 public:
  /**
   * @brief Construct a new Node object
   *
   * @param id Global ID of the node
   * @param inlet_eles Inlet element of the node
   * @param outlet_eles Outlet element of the node
   * @param model The model to which the node belongs
   */
  Node(int id, const std::vector<Block *> &inlet_eles,
       const std::vector<Block *> &outlet_eles, Model *model);

  int id;                            ///< Global ID of the block
  std::vector<Block *> inlet_eles;   ///< Inlet element of the node
  std::vector<Block *> outlet_eles;  ///< Outlet element of the node
  Model *model{nullptr};             ///< The model to which the node belongs

  int flow_dof{0};  ///< Global flow degree-of-freedom of the node
  int pres_dof{0};  ///< Global pressure degree-of-freedom of the node

  /**
   * @brief Get the name of the node
   *
   * @return std::string Name of the node
   */
  std::string get_name();

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set global_var_ids and global_eqn_ids of the element based on the
   * number of equations and the number of internal variables of the
   * element.
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler &dofhandler);
};

#endif  // SVZERODSOLVER_MODEL_NODE_HPP_
