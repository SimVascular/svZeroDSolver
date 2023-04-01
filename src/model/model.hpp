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
 * @file model.hpp
 * @brief MODEL::Model source file
 */
#ifndef SVZERODSOLVER_MODEL_MODEL_HPP_
#define SVZERODSOLVER_MODEL_MODEL_HPP_

#include <algorithm>
#include <list>
#include <memory>
#include <vector>

#include "../algebra/sparsesystem.hpp"
#include "../model/bloodvessel.hpp"
#include "../model/bloodvesseljunction.hpp"
#include "../model/closedloopRCRbc.hpp"
#include "../model/closedloopcoronarybc.hpp"
#include "../model/closedloopheartpulmonary.hpp"
#include "../model/flowreferencebc.hpp"
#include "../model/junction.hpp"
#include "../model/openloopcoronarybc.hpp"
#include "../model/pressurereferencebc.hpp"
#include "../model/resistancebc.hpp"
#include "../model/resistivejunction.hpp"
#include "../model/windkesselbc.hpp"
#include "block.hpp"
#include "blocktype.hpp"
#include "dofhandler.hpp"
#include "node.hpp"
#include "parameter.hpp"

namespace MODEL {

/**
 * @brief Model of 0D elements
 *
 * This class represents a full 0D model. It contains attributes and
 * methods to store and modify 0D elements.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class Model {
 public:
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

  T cardiac_cycle_period = -1.0;  ///< Cardiac cycle period
  T time = 0.0;                   ///< Current time

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
                std::string_view name, bool internal = false);

  /**
   * @brief Get a block by its name
   *
   * @param name Name of the Block
   * @return Block<T>* The block
   */
  Block<T> *get_block(std::string_view name);

  /**
   * @brief Get a block by its global ID
   *
   * @param block_id Global ID of the Block
   * @return Block<T>* The block
   */
  Block<T> *get_block(int block_id);

  /**
   * @brief Get the name of a block by it's ID
   *
   * @param block_id Global ID of the block
   * @return std::string Name of the block
   */
  std::string get_block_name(int block_id);

  /**
   * @brief Add a node to the model
   *
   * @param inlet_eles Inlet blocks of the node
   * @param outlet_eles Outlet blocks of the node
   * @return int Global ID of the node
   */
  int add_node(const std::vector<Block<T> *> &inlet_eles,
               const std::vector<Block<T> *> &outlet_eles,
               std::string_view name);

  /**
   * @brief Get the name of a node by it's ID
   *
   * @param node_id Global ID of the node
   * @return std::string Name of the node
   */
  std::string get_node_name(int node_id);

  /**
   * @brief Add a constant model parameter
   *
   * @param value Value of the parameter
   * @return int Global ID of the parameter
   */
  int add_parameter(T value);

  /**
   * @brief Add a time-dependent model parameter
   *
   * @param times Times corresponding to the parameter values
   * @param values Values of the parameter
   * @param periodic Toggle whether parameter is periodic
   * @return int Global ID of the parameter
   */
  int add_parameter(const std::vector<T> &times, const std::vector<T> &values,
                    bool periodic = true);

  /**
   * @brief Get a parameter by its global ID
   *
   * @param param_id Global ID of the parameter
   * @return Parameter<T>* The parameter
   */
  Parameter<T> *get_parameter(int param_id);

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
  void update_constant(ALGEBRA::SparseSystem<T> &system);

  /**
   * @brief Update the time-dependent contributions of all elements in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param time Current time
   */
  void update_time(ALGEBRA::SparseSystem<T> &system, T time);

  /**
   * @brief Update the solution-dependent contributions of all elements in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(ALGEBRA::SparseSystem<T> &system,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &dy);

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
  std::map<std::string, int> get_num_triplets();

  /**
   * @brief Get the number of blocks in the model
   *
   * @return int Number of blocks
   */
  int get_num_blocks(bool internal = false);

 private:
  int block_count = 0;
  int node_count = 0;
  int parameter_count = 0;
  std::map<int, T> param_value_cache;

  std::vector<std::shared_ptr<Block<T>>> blocks;  ///< Blocks of the model
  std::vector<BlockType> block_types;             ///< Types of the blocks
  std::vector<std::string> block_names;           ///< Names of the blocks
  std::map<std::string_view, int>
      block_index_map;  ///< Map between block name and index

  std::vector<std::shared_ptr<Block<T>>>
      hidden_blocks;  ///< Hidden blocks of the model

  std::vector<std::shared_ptr<Node<T>>> nodes;  ///< Nodes of the model
  std::vector<std::string> node_names;          ///< Names of the nodes

  std::vector<std::shared_ptr<Parameter<T>>>
      parameters;  ///< Parameters of the model
  std::vector<std::shared_ptr<Parameter<T>>>
      time_dependent_parameters;  ///< Time-dependent parameters of the model
  std::vector<double> parameter_values;  ///< Current values of the parameters
};

template <typename T>
Model<T>::Model() {}

template <typename T>
Model<T>::~Model() {}

template <typename T>
int Model<T>::add_block(BlockType block_type,
                        const std::vector<int> &block_param_ids,
                        std::string_view name, bool internal) {
  DEBUG_MSG("Adding block " << name << " with type " << block_type);
  std::shared_ptr<Block<T>> block;
  switch (block_type) {
    case BlockType::BLOODVESSEL:
      block = std::shared_ptr<Block<T>>(
          new BloodVessel<T>(block_count, block_param_ids, this));
      break;
    case BlockType::JUNCTION:
      block = std::shared_ptr<Block<T>>(
          new Junction<T>(block_count, block_param_ids, this));
      break;
    case BlockType::BLOODVESSELJUNCTION:
      block = std::shared_ptr<Block<T>>(
          new BloodVesselJunction<T>(block_count, block_param_ids, this));
      break;
    case BlockType::RESISTIVEJUNCTION:
      block = std::shared_ptr<Block<T>>(
          new ResistiveJunction<T>(block_count, block_param_ids, this));
      break;
    case BlockType::FLOWBC:
      block = std::shared_ptr<Block<T>>(
          new FlowReferenceBC<T>(block_count, block_param_ids, this));
      break;
    case BlockType::RESISTANCEBC:
      block = std::shared_ptr<Block<T>>(
          new ResistanceBC<T>(block_count, block_param_ids, this));
      break;
    case BlockType::WINDKESSELBC:
      block = std::shared_ptr<Block<T>>(
          new WindkesselBC<T>(block_count, block_param_ids, this));
      break;
    case BlockType::PRESSUREBC:
      block = std::shared_ptr<Block<T>>(
          new PressureReferenceBC<T>(block_count, block_param_ids, this));
      break;
    case BlockType::OPENLOOPCORONARYBC:
      block = std::shared_ptr<Block<T>>(
          new OpenLoopCoronaryBC<T>(block_count, block_param_ids, this));
      break;
    case BlockType::CLOSEDLOOPCORONARYLEFTBC:
      block = std::shared_ptr<Block<T>>(
          new ClosedLoopCoronaryBC<T, MODEL::Side::LEFT>(
              block_count, block_param_ids, this));
      break;
    case BlockType::CLOSEDLOOPCORONARYRIGHTBC:
      block = std::shared_ptr<Block<T>>(
          new ClosedLoopCoronaryBC<T, MODEL::Side::RIGHT>(
              block_count, block_param_ids, this));
      break;
    case BlockType::CLOSEDLOOPRCRBC:
      block = std::shared_ptr<Block<T>>(
          new ClosedLoopRCRBC<T>(block_count, block_param_ids, this));
      break;
    case BlockType::CLOSEDLOOPHEARTPULMONARY:
      block = std::shared_ptr<Block<T>>(
          new ClosedLoopHeartPulmonary<T>(block_count, block_param_ids, this));
      break;
    default:
      throw std::runtime_error(
          "Adding block to model failed: Invalid block type!");
  }
  if (internal) {
    hidden_blocks.push_back(block);
  } else {
    blocks.push_back(block);
    block_types.push_back(block_type);
    block_index_map.insert({name, block_count});
    block_names.push_back(static_cast<std::string>(name));
  }
  return block_count++;
}

template <typename T>
Block<T> *Model<T>::get_block(std::string_view name) {
  return blocks[block_index_map[name]].get();
}

template <typename T>
Block<T> *Model<T>::get_block(int block_id) {
  if (block_id >= blocks.size()) {
    return hidden_blocks[block_id - blocks.size()].get();
  }
  return blocks[block_id].get();
}

template <typename T>
std::string Model<T>::get_block_name(int block_id) {
  return block_names[block_id];
}

template <typename T>
int Model<T>::add_node(const std::vector<Block<T> *> &inlet_eles,
                       const std::vector<Block<T> *> &outlet_eles,
                       std::string_view name) {
  DEBUG_MSG("Adding node " << name);
  auto node = std::shared_ptr<Node<T>>(
      new Node<T>(node_count, inlet_eles, outlet_eles, this));
  nodes.push_back(node);
  node_names.push_back(static_cast<std::string>(name));
  return node_count++;
}

template <typename T>
std::string Model<T>::get_node_name(int node_id) {
  return node_names[node_id];
}

template <typename T>
int Model<T>::add_parameter(T value) {
  std::shared_ptr<Parameter<T>> param(new Parameter<T>(parameter_count, value));
  parameters.push_back(param);
  parameter_values.push_back(param->get(0.0));
  return parameter_count++;
}

template <typename T>
int Model<T>::add_parameter(const std::vector<T> &times,
                            const std::vector<T> &values, bool periodic) {
  std::shared_ptr<Parameter<T>> param(
      new Parameter<T>(parameter_count, times, values, periodic));
  if (periodic && (param->isconstant == false)) {
    if ((this->cardiac_cycle_period > 0.0) &&
        (param->cycle_period != this->cardiac_cycle_period)) {
      throw std::runtime_error(
          "Inconsistent cardiac cycle period defined in parameters");
    }
    this->cardiac_cycle_period = param->cycle_period;
  }
  parameters.push_back(param);
  time_dependent_parameters.push_back(param);
  parameter_values.push_back(param->get(0.0));
  return parameter_count++;
}

template <typename T>
Parameter<T> *Model<T>::get_parameter(int param_id) {
  return parameters[param_id].get();
}

template <typename T>
void Model<T>::finalize() {
  for (auto &node : nodes) {
    node->setup_dofs(dofhandler);
  }
  for (auto &block : blocks) {
    block->setup_dofs(dofhandler);
  }
  for (auto &block : blocks) {
    block->setup_model_dependent_params();
  }
}

template <typename T>
int Model<T>::get_num_blocks(bool internal) {
  int num_blocks = blocks.size();
  if (internal) {
    num_blocks += hidden_blocks.size();
  }
  return num_blocks;
}

template <typename T>
void Model<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  for (auto block : blocks) {
    block->update_constant(system, parameter_values);
  }
}

template <typename T>
void Model<T>::update_time(ALGEBRA::SparseSystem<T> &system, T time) {
  this->time = time;
  for (auto param : time_dependent_parameters) {
    parameter_values[param->id] = param->get(time);
  }
  for (auto block : blocks) {
    block->update_time(system, parameter_values);
  }
}

template <typename T>
void Model<T>::update_solution(ALGEBRA::SparseSystem<T> &system,
                               Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
                               Eigen::Matrix<T, Eigen::Dynamic, 1> &dy) {
  for (auto block : blocks) {
    block->update_solution(system, parameter_values, y, dy);
  }
}

template <typename T>
void Model<T>::to_steady() {
  for (auto param : time_dependent_parameters) {
    param->to_steady();
  }
  for (size_t i = 0; i < block_types.size(); i++) {
    blocks[i]->steady = true;
    if ((block_types[i] == BlockType::WINDKESSELBC) ||
        (block_types[i] == BlockType::CLOSEDLOOPRCRBC)) {
      int param_id_capacitance = blocks[i]->global_param_ids[1];
      T value = parameters[param_id_capacitance]->get(0.0);
      param_value_cache.insert({param_id_capacitance, value});
    }
  }
}

template <typename T>
void Model<T>::to_unsteady() {
  for (auto param : time_dependent_parameters) {
    param->to_unsteady();
  }
  for (auto &[param_id_capacitance, value] : param_value_cache) {
    DEBUG_MSG("Setting Windkessel capacitance back to " << value);
    parameters[param_id_capacitance]->update(value);
  }
  for (size_t i = 0; i < block_types.size(); i++) {
    blocks[i]->steady = false;
  }
}

template <typename T>
std::map<std::string, int> Model<T>::get_num_triplets() {
  std::map<std::string, int> num_triplets = {
      {"F", 0},
      {"E", 0},
      {"D", 0},
  };
  for (auto &elem : blocks) {
    for (auto &[key, value] : elem->get_num_triplets()) {
      num_triplets[key] += value;
    }
  }
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_MODEL_HPP_
