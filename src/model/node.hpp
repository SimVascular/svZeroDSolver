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
  flow_dof = dofhandler.register_variable("Q_" + name);
  pres_dof = dofhandler.register_variable("P_" + name);
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_NODE_HPP_