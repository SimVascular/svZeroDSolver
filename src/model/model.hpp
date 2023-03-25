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
#include "block.hpp"
#include "dofhandler.hpp"
#include "node.hpp"
#include "parameter.hpp"
// #include "../model/closedloopRCRbc.hpp"
// #include "../model/closedloopcoronarybc.hpp"
// #include "../model/closedloopheartpulmonary.hpp"
#include "../model/flowreferencebc.hpp"
#include "../model/junction.hpp"
#include "../model/windkesselbc.hpp"
// #include "../model/openloopcoronarybc.hpp"
#include "../model/pressurereferencebc.hpp"
#include "../model/resistancebc.hpp"
// #include "../model/resistivejunction.hpp

namespace MODEL {

enum BlockType {
  BLOODVESSEL = 0,
  JUNCTION = 1,
  BLOODVESSELJUNCTION = 2,
  FLOWBC = 3,
  PRESSUREBC = 4,
  RESISTANCEBC = 5,
  WINDKESSELBC = 6,
};

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

  std::vector<std::shared_ptr<Block<T>>> blocks;  ///< Blocks of the model
  std::vector<BlockType> block_types;             ///< Types of the blocks
  std::vector<std::shared_ptr<Parameter<T>>>
      parameters;  ///< Parameters of the model
  std::vector<std::shared_ptr<Parameter<T>>>
      time_dependent_parameters;  ///< Time-dependent parameters of the model
  std::vector<double> parameter_values;  ///< Current values of the parameters
  std::map<std::string_view, int>
      block_index_map;        ///< Map between block name and index
  DOFHandler dofhandler;      ///< Degree-of-freedom handler of the model
  std::vector<Node *> nodes;  ///< Nodes of the model
  std::vector<std::string>
      external_coupling_blocks;  ///< List of external coupling blocks (names
                                 ///< need to be available for interface)

  /**
   * @brief Add a block to the model
   *
   * @param block_type Type of the block
   * @param block_param_ids Global IDs of the parameters of the block
   * @param name The name of the block
   */
  int add_block(BlockType block_type, const std::vector<int> &block_param_ids,
                std::string_view name);

  /**
   * @brief Add a constant model parameter
   *
   * @param value Value of the parameter
   */
  int add_parameter(T value);

  /**
   * @brief Add a time-dependent model parameter
   *
   * @param times Times corresponding to the parameter values
   * @param values Values of the parameter
   * @param periodic Toggle whether parameter is periodic
   */
  int add_parameter(const std::vector<T> &times, const std::vector<T> &values,
                    bool periodic = true);

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
   */
  void update_solution(ALGEBRA::SparseSystem<T> &system,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

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

 private:
  int block_count = 0;
  int parameter_count = 0;
  std::map<int, T> param_value_cache;
};

template <typename T>
Model<T>::Model() {}

template <typename T>
Model<T>::~Model() {}

template <typename T>
int Model<T>::add_block(BlockType block_type,
                        const std::vector<int> &block_param_ids,
                        std::string_view name) {
  DEBUG_MSG("Adding block " << name << " with type " << block_type);
  std::shared_ptr<Block<T>> block;
  switch (block_type) {
    case BlockType::BLOODVESSEL:
      block = std::shared_ptr<Block<T>>(
          new BloodVessel<T>(block_count, block_param_ids, name));
      break;
    case BlockType::JUNCTION:
      block = std::shared_ptr<Block<T>>(
          new Junction<T>(block_count, block_param_ids, name));
      break;
    case BlockType::FLOWBC:
      block = std::shared_ptr<Block<T>>(
          new FlowReferenceBC<T>(block_count, block_param_ids, name));
      break;
    case BlockType::RESISTANCEBC:
      block = std::shared_ptr<Block<T>>(
          new ResistanceBC<T>(block_count, block_param_ids, name));
      break;
    case BlockType::WINDKESSELBC:
      block = std::shared_ptr<Block<T>>(
          new WindkesselBC<T>(block_count, block_param_ids, name));
      break;
    case BlockType::PRESSUREBC:
      block = std::shared_ptr<Block<T>>(
          new PressureReferenceBC<T>(block_count, block_param_ids, name));
      break;
    case BlockType::BLOODVESSELJUNCTION:
      block = std::shared_ptr<Block<T>>(
          new BloodVesselJunction<T>(block_count, block_param_ids, name));
      break;
    default:
      throw std::runtime_error(
          "Adding block to model failed: Invalid block type!");
  }
  blocks.push_back(block);
  block_types.push_back(block_type);
  block_index_map.insert({name, block_count});
  return block_count++;
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
  parameters.push_back(param);
  time_dependent_parameters.push_back(param);
  parameter_values.push_back(param->get(0.0));
  return parameter_count++;
}

template <typename T>
void Model<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  for (auto block : blocks) {
    block->update_constant(system, parameter_values);
  }
}

template <typename T>
void Model<T>::update_time(ALGEBRA::SparseSystem<T> &system, T time) {
  for (auto param : time_dependent_parameters) {
    parameter_values[param->id] = param->get(time);
  }
  for (auto block : blocks) {
    block->update_time(system, parameter_values);
  }
}

template <typename T>
void Model<T>::update_solution(ALGEBRA::SparseSystem<T> &system,
                               Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  for (auto block : blocks) {
    block->update_solution(system, parameter_values, y);
  }
}

template <typename T>
void Model<T>::to_steady() {
  for (auto param : time_dependent_parameters) {
    param->to_steady();
  }
  for (size_t i = 0; i < block_types.size(); i++) {
    if (block_types[i] == BlockType::WINDKESSELBC) {
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
