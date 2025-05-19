// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
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
