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
 * @file Model.h
 * @brief model::Model source file
 */

#ifndef SVZERODSOLVER_MODEL_MODEL_HPP_
#define SVZERODSOLVER_MODEL_MODEL_HPP_

#include <algorithm>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Block.h"
#include "BlockFactory.h"
#include "BloodVessel.h"
#include "BloodVesselJunction.h"
#include "ClosedLoopCoronaryBC.h"
#include "ClosedLoopHeartPulmonary.h"
#include "ClosedLoopRCRBC.h"
#include "DOFHandler.h"
#include "FlowReferenceBC.h"
#include "Junction.h"
#include "Node.h"
#include "OpenLoopCoronaryBC.h"
#include "Parameter.h"
#include "PressureReferenceBC.h"
#include "ResistanceBC.h"
#include "ResistiveJunction.h"
#include "WindkesselBC.h"

/**
 * @brief Model of 0D elements
 *
 * This class represents a full 0D model. It contains attributes and
 * methods to store and modify 0D elements.
 *
 */
class Model {
 public:
  /// Factory that holds all implemented blocks
  std::map<BlockType, BlockFactoryFunc> blockFactoryMap;

  /**
   * @brief Construct a new Model object
   *
   */
  Model();

  /**
   * @brief Destroy the Model object
   *
   */
  ~Model();

  DOFHandler dofhandler;  ///< Degree-of-freedom handler of the model

  double cardiac_cycle_period = -1.0;  ///< Cardiac cycle period
  double time = 0.0;                   ///< Current time

  /**
   * @brief Add a block to the model
   *
   * @param block_type Type of the block
   * @param block_param_ids Global IDs of the parameters of the block
   * @param name The name of the block
   * @param internal Toggle whether block is internal
   * @return int Global ID of the block
   */
  int add_block(BlockType block_type, const std::vector<int> &block_param_ids,
                const std::string_view &name, bool internal = false);

  /**
   * @brief Get a block by its name
   *
   * @param name Name of the Block
   * @return Block* The block
   */
  Block *get_block(const std::string_view &name) const;

  /**
   * @brief Get a block by its global ID
   *
   * @param block_id Global ID of the Block
   * @return Block* The block
   */
  Block *get_block(int block_id) const;

  /**
   * @brief Get a block type by its name
   *
   * @param block_id Global ID of the Block
   * @return BlockType The block type
   */
  BlockType get_block_type(const std::string_view &name) const;

  /**
   * @brief Get the name of a block by it's ID
   *
   * @param block_id Global ID of the block
   * @return std::string Name of the block
   */
  std::string get_block_name(int block_id) const;

  /**
   * @brief Add a node to the model
   *
   * @param inlet_eles Inlet blocks of the node
   * @param outlet_eles Outlet blocks of the node
   * @return int Global ID of the node
   */
  int add_node(const std::vector<Block *> &inlet_eles,
               const std::vector<Block *> &outlet_eles,
               const std::string_view &name);

  /**
   * @brief Get the name of a node by it's ID
   *
   * @param node_id Global ID of the node
   * @return std::string Name of the node
   */
  std::string get_node_name(int node_id) const;

  /**
   * @brief Add a constant model parameter
   *
   * @param value Value of the parameter
   * @return int Global ID of the parameter
   */
  int add_parameter(double value);

  /**
   * @brief Add a time-dependent model parameter
   *
   * @param times Times corresponding to the parameter values
   * @param values Values of the parameter
   * @param periodic Toggle whether parameter is periodic
   * @return int Global ID of the parameter
   */
  int add_parameter(const std::vector<double> &times,
                    const std::vector<double> &values, bool periodic = true);

  /**
   * @brief Get a parameter by its global ID
   *
   * @param param_id Global ID of the parameter
   * @return Parameter* The parameter
   */
  Parameter *get_parameter(int param_id);

  /**
   * @brief Get the current value of a parameter
   *
   * @param param_id Global ID of the parameter
   * @return T Current value of the parameter
   */
  double get_parameter_value(int param_id) const;

  /**
   * @brief Update the current value of a parameter
   *
   * @param param_id Global ID of the parameter
   * @param param_value The new parameter value
   */
  void update_parameter_value(int param_id, double param_value);

  /**
   * @brief Finalize the model after all blocks, nodes and parameters have been
   * added
   *
   */
  void finalize();

  /**
   * @brief Update the constant contributions of all elements in a sparse system
   *
   * @param system System to update contributions at
   */
  void update_constant(SparseSystem &system);

  /**
   * @brief Update the time-dependent contributions of all elements in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param time Current time
   */
  void update_time(SparseSystem &system, double time);

  /**
   * @brief Update the solution-dependent contributions of all elements in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(SparseSystem &system,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  /**
   * @brief Convert the blocks to a steady behavior
   *
   */
  void to_steady();

  /**
   * @brief Convert the blocks to an unsteady behavior
   *
   */
  void to_unsteady();

  /**
   * @brief Get number of triplets all elements
   *
   * Get the number of triplets the elements contribute to the global system
   * (relevant for sparse memory reservation)
   *
   * @return Number of triplets that are used in each system matrix
   */
  // std::map<std::string, int> get_num_triplets();
  TripletsContributions get_num_triplets() const;

  /**
   * @brief Get the number of blocks in the model
   *
   * @return int Number of blocks
   */
  int get_num_blocks(bool internal = false) const;

 private:
  int block_count = 0;
  int node_count = 0;
  int parameter_count = 0;
  std::map<int, double> param_value_cache;

  std::vector<std::shared_ptr<Block>> blocks;  ///< Blocks of the model
  std::vector<BlockType> block_types;          ///< Types of the blocks
  std::vector<std::string> block_names;        ///< Names of the blocks
  std::map<std::string, int>
      block_index_map;  ///< Map between block name and index

  std::vector<std::shared_ptr<Block>>
      hidden_blocks;  ///< Hidden blocks of the model

  std::vector<std::shared_ptr<Node>> nodes;  ///< Nodes of the model
  std::vector<std::string> node_names;       ///< Names of the nodes

  std::vector<Parameter> parameters;     ///< Parameters of the model
  std::vector<double> parameter_values;  ///< Current values of the parameters
};

#endif  // SVZERODSOLVER_MODEL_MODEL_HPP_
