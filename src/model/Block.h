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
 * @file Block.h
 * @brief model::Block source file
 */
#ifndef SVZERODSOLVER_MODEL_BLOCK_HPP_
#define SVZERODSOLVER_MODEL_BLOCK_HPP_

#include <Eigen/Core>
#include <list>
#include <map>
#include <vector>

#include "DOFHandler.h"
#include "SparseSystem.h"

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
   * @brief Contributions to D matrix
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
  /**
   * @brief Construct a new Block object
   *
   * @param id Global ID of the block
   * @param param_ids Global IDs of the block parameters
   * @param model The model to which the block belongs
   */
  explicit Block(int id, const std::vector<int> &param_ids, Model *model);

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

  int id;        ///< Global ID of the block
  Model *model;  ///< The model to which the block belongs

  std::vector<Node *> inlet_nodes;   ///< Inlet nodes
  std::vector<Node *> outlet_nodes;  ///< Outlet nodes

  bool steady = false;  ///< Toggle steady behavior

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
  virtual void update_solution(SparseSystem &system,
                               std::vector<double> &parameters,
                               Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                               Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

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
