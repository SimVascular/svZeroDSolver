/**
 * @file block.hpp
 * @brief MODEL::Block source file
 */
#ifndef SVZERODSOLVER_MODEL_BLOCK_HPP_
#define SVZERODSOLVER_MODEL_BLOCK_HPP_

#include <map>
#include <vector>

#include "../algebra/densesystem.hpp"
#include "../algebra/sparsesystem.hpp"
#include "dofhandler.hpp"
#include "node.hpp"

namespace MODEL {

/**
 * @brief Base class for 0D model components.
 *
 * A MODEL::Block is the base class of 0D model elements. It is the place
 * where the contribution of an element to the global system is controlled.
 * It defines all relevant attributes and methods of an element and a few
 * common helpers for setting it up.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class Block {
 public:
  /**
   * @brief Parameters of the element.
   *
   * Struct containing all constant and/or time-dependent parameters of the
   * element.
   */
  struct Parameters {};

  /**
   * @brief Construct a new Block object
   *
   */
  Block();

  /**
   * @brief Construct a new Block object
   *
   * @param name Name of the block
   */
  Block(std::string name);

  /**
   * @brief Destroy the Block object
   *
   */
  ~Block();

  std::string name;                  ///< Name of the block
  std::vector<Node *> inlet_nodes;   ///< Inlet nodes
  std::vector<Node *> outlet_nodes;  ///< Outlet nodes

  /**
   * @brief Global variable indices of the local element contributions
   *
   * Determines where the local element contributions are written to in
   * the global system during assembly. The order of indices is:
   *
   * \f[
   * [P_{in}^1, Q_{in}^1, \dots, P_{in}^n, Q_{in}^n, P_{out}^1, Q_{out}^1,
   * \dots, P_{out}^m, Q_{out}^m, V^{1}, \dots, V^{p}] \f]
   *
   * with \f$P_{in} \f$, \f$Q_{in} \f$, \f$P_{out} \f$, \f$Q_{out} \f$, and \f$V
   * \f$ denoting inlet pressure, inlet flow, outlet pressure, outlet flow and
   * an internal variable of the element, respectively.
   *
   * Variable indices correspond to columns in the global system.
   *
   */
  std::vector<unsigned int> global_var_ids;

  /**
   * @brief Global equation indices of the local element contributions
   *
   * Equation indices correspond to rows in the global system.
   */
  std::vector<unsigned int> global_eqn_ids;

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set \ref global_var_ids and \ref global_eqn_ids of the element based on the
   * number of equations and the number of internal variables of the
   * element.
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   * @param num_equations Number of equations of the block
   * @param num_internal_vars Number of internal variables of the block
   */
  void setup_dofs_(DOFHandler &dofhandler, unsigned int num_equations,
                   unsigned int num_internal_vars);

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

  /**
   * @brief Update the constant contributions of the element in a dense system
   *
   * @param system System to update contributions at
   */
  void update_constant(ALGEBRA::DenseSystem<T> &system);

  /**
   * @brief Update the time-dependent contributions of the element in a dense
   * system
   *
   * @param system System to update contributions at
   * @param time Current time
   */
  void update_time(ALGEBRA::DenseSystem<T> &system, T time);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * dense system
   *
   * @param system System to update contributions at
   * @param y Current solution
   */
  void update_solution(ALGEBRA::DenseSystem<T> &system,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   */
  void update_constant(ALGEBRA::SparseSystem<T> &system);

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param time Current time
   */
  void update_time(ALGEBRA::SparseSystem<T> &system, T time);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param y Current solution
   */
  void update_solution(ALGEBRA::SparseSystem<T> &system,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

  /**
   * @brief Convert the block to a steady behavior
   *
   */
  void to_steady();

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  std::map<std::string, int> num_triplets = {
      {"F", 0},
      {"E", 0},
      {"D", 0},
  };

 private:
  Parameters params;  ///< Parameters of the element
};

template <typename T>
Block<T>::Block() {}

template <typename T>
Block<T>::Block(std::string name) {
  this->name = name;
}

template <typename T>
Block<T>::~Block() {}

template <typename T>
void Block<T>::setup_dofs_(DOFHandler &dofhandler, unsigned int num_equations,
                           unsigned int num_internal_vars) {
  // Collect external DOFs from inlet nodes
  for (auto inlet_node : inlet_nodes) {
    global_var_ids.push_back(inlet_node->pres_dof);
    global_var_ids.push_back(inlet_node->flow_dof);
  }

  // Collect external DOFs from outlet nodes
  for (auto outlet_node : outlet_nodes) {
    global_var_ids.push_back(outlet_node->pres_dof);
    global_var_ids.push_back(outlet_node->flow_dof);
  }

  // Register internal variables of block
  for (unsigned int i = 0; i < num_internal_vars; i++) {
    global_var_ids.push_back(
        dofhandler.register_variable("var_" + std::to_string(i) + "_" + name));
  }

  // Register equations of block
  for (unsigned int i = 0; i < num_equations; i++) {
    global_eqn_ids.push_back(dofhandler.register_equation());
  }
}

template <typename T>
void Block<T>::setup_dofs(DOFHandler &dofhandler) {}

template <typename T>
void Block<T>::update_constant(ALGEBRA::DenseSystem<T> &system) {}

template <typename T>
void Block<T>::update_time(ALGEBRA::DenseSystem<T> &system, T time) {}

template <typename T>
void Block<T>::update_solution(ALGEBRA::DenseSystem<T> &system,
                               Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {}

template <typename T>
void Block<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {}

template <typename T>
void Block<T>::update_time(ALGEBRA::SparseSystem<T> &system, T time) {}

template <typename T>
void Block<T>::update_solution(ALGEBRA::SparseSystem<T> &system,
                               Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {}

template <typename T>
void Block<T>::to_steady() {}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_BLOCK_HPP_