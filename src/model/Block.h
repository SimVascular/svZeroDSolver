// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file Block.h
 * @brief model::Block source file
 */
#ifndef SVZERODSOLVER_MODEL_BLOCK_HPP_
#define SVZERODSOLVER_MODEL_BLOCK_HPP_

#include <Eigen/Core>
#include <list>
#include <map>
#include <vector>

#include "BlockType.h"
#include "DOFHandler.h"
#include "Parameter.h"
#include "SparseSystem.h"
#include "State.h"

/**
 * @brief The number of triplets that the element contributes
 * to the global system.
 */
struct TripletsContributions {
  TripletsContributions(){};
  /**
   * @brief Set the number of triplets that the element contributes
   * to the global system.
   * @param F Contributions to F matrix
   * @param E Contributions to E matrix
   * @param D Contributions to dC/dy matrix
   */
  TripletsContributions(int F, int E, int D) : F(F), E(E), D(D){};
  /**
   * @brief Set the number of triplets that the element contributes
   * to the global system.
   * @param other TripletsContributions object to add to the
   * number of contributions
   * @return The number of triplets
   */
  TripletsContributions operator+=(const TripletsContributions &other) {
    F += other.F;
    E += other.E;
    D += other.D;
    return *this;
  };

  /**
   * @brief Contributions to F matrix
   */
  int F{0};
  /**
   * @brief Contributions to E matrix
   */
  int E{0};
  /**
   * @brief Contributions to dC/dy matrix
   */
  int D{0};
};

class Node;
class Model;

/**
 * @brief Base class for 0D model components.
 *
 * A Block is the base class of 0D model elements. It is the place
 * where the contribution of an element to the global system is controlled.
 * It defines all relevant attributes and methods of an element and a few
 * common helpers for setting it up.
 */
class Block {
 public:
  const int id;                  ///< Global ID of the block
  const Model *model;            ///< The model to which the block belongs
  const BlockType block_type;    ///< Type of this block
  const BlockClass block_class;  ///< Class of this block
  VesselType vessel_type = VesselType::neither;  ///< Vessel type of this block
  const std::vector<std::pair<std::string, InputParameter>>
      input_params;  ///< Map from name to input parameter

  std::vector<Node *> inlet_nodes;   ///< Inlet nodes
  std::vector<Node *> outlet_nodes;  ///< Outlet nodes

  bool steady = false;             ///< Toggle steady behavior
  bool input_params_list = false;  ///< Are input parameters given as a list?

  /**
   * @brief Construct a new Block object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   * @param block_type The specific type of block
   * @param block_class The class the block belongs to (e.g. vessel, junction)
   * @param input_params The parameters the block takes from the input file
   */
  Block(int id, Model *model, BlockType block_type, BlockClass block_class,
        std::vector<std::pair<std::string, InputParameter>> input_params)
      : id(id),
        model(model),
        block_type(block_type),
        block_class(block_class),
        input_params(input_params) {}

  /**
   * @brief Destroy the Block object
   *
   */
  ~Block();

  /**
   * @brief Copy the Block object
   *
   */
  Block(const Block &) = delete;
  /**
   * @brief Global IDs for the block parameters.
   *
   */
  std::vector<int> global_param_ids;  ///< IDs of the parameters

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
  std::vector<int> global_var_ids;

  /**
   * @brief Global equation indices of the local element contributions
   *
   * Equation indices correspond to rows in the global system.
   */
  std::vector<int> global_eqn_ids;

  /**
   * @brief Get the name of the block
   *
   * @return std::string Name of the block
   */
  std::string get_name();

  /**
   * @brief Update vessel type of the block
   *
   * @param type Type of vessel
   */
  void update_vessel_type(VesselType type);

  /**
   * @brief Setup parameter IDs for the block
   * @param param_ids Global IDs of the block parameters
   */
  void setup_params_(const std::vector<int> &param_ids);

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
   * @param internal_var_names Number of internal variables of the block
   */

  void setup_dofs_(DOFHandler &dofhandler, int num_equations,
                   const std::list<std::string> &internal_var_names);

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
  virtual void setup_dofs(DOFHandler &dofhandler);

  /**
   * @brief Setup parameters that depend on the model
   *
   */
  virtual void setup_model_dependent_params();

  /**
   * @brief Setup parameters that depend on the initial state
   *
   * @param initial_state The initial state of the system
   * @param parameters The parameter values vector (at time 0)
   */
  virtual void setup_initial_state_dependent_params(
      State initial_state, std::vector<double> &parameters);

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  virtual void update_constant(SparseSystem &system,
                               std::vector<double> &parameters);
  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  virtual void update_time(SparseSystem &system,
                           std::vector<double> &parameters);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  virtual void update_solution(
      SparseSystem &system, std::vector<double> &parameters,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  /**
   * @brief Modify the solution after solving it
   *
   * @param y Current solution
   */
  virtual void post_solve(Eigen::Matrix<double, Eigen::Dynamic, 1> &y);

  /**
   * @brief Set the gradient of the block contributions with respect to the
   * parameters
   *
   * @param jacobian Jacobian with respect to the parameters
   * @param alpha Current parameter vector
   * @param residual Residual with respect to the parameters
   * @param y Current solution
   * @param dy Time-derivative of the current solution
   */
  virtual void update_gradient(
      Eigen::SparseMatrix<double> &jacobian,
      Eigen::Matrix<double, Eigen::Dynamic, 1> &residual,
      Eigen::Matrix<double, Eigen::Dynamic, 1> &alpha, std::vector<double> &y,
      std::vector<double> &dy);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets;

  /**
   * @brief Get number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   *
   * @return TripletsContributions Number of triplets of element
   */
  virtual TripletsContributions get_num_triplets();
};

#endif
