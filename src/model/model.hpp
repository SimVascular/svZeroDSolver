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

#include <list>
#include <map>
#include <variant>

#include "../algebra/densesystem.hpp"
#include "../algebra/sparsesystem.hpp"
#include "bloodvessel.hpp"
#include "dofhandler.hpp"
#include "flowreferencebc.hpp"
#include "junction.hpp"
#include "node.hpp"
#include "openloopcoronarybc.hpp"
#include "pressurereferencebc.hpp"
#include "resistancebc.hpp"
#include "windkesselbc.hpp"

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

  std::map<std::string_view,
           std::variant<Junction<T>, BloodVessel<T>, FlowReferenceBC<T>,
                        PressureReferenceBC<T>, WindkesselBC<T>,
                        ResistanceBC<T>, OpenLoopCoronaryBC<T>>>
      blocks;             ///< Elements of the model
  DOFHandler dofhandler;  ///< Degree-of-freedom handler of the model
  std::list<Node> nodes;  ///< Nodes of the model

  /**
   * @brief Update the constant contributions of all elements in a dense system
   *
   * @param system System to update contributions at
   */
  void update_constant(ALGEBRA::DenseSystem<T> &system);

  /**
   * @brief Update the time-dependent contributions of all elements in a dense
   * system
   *
   * @param system System to update contributions at
   * @param time Current time
   */
  void update_time(ALGEBRA::DenseSystem<T> &system, T time);

  /**
   * @brief Update the solution-dependent contributions of all elements in a
   * dense system
   *
   * @param system System to update contributions at
   * @param y Current solution
   */
  void update_solution(ALGEBRA::DenseSystem<T> &system,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

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
};

template <typename T>
Model<T>::Model() {}

template <typename T>
Model<T>::~Model() {}

template <typename T>
void Model<T>::update_constant(ALGEBRA::DenseSystem<T> &system) {
  for (auto &&elem : blocks) {
    std::visit([&](auto &&block) { block.update_constant(system); },
               elem.second);
  }
}

template <typename T>
void Model<T>::update_time(ALGEBRA::DenseSystem<T> &system, T time) {
  for (auto &&elem : blocks) {
    std::visit([&](auto &&block) { block.update_time(system, time); },
               elem.second);
  }
}

template <typename T>
void Model<T>::update_solution(ALGEBRA::DenseSystem<T> &system,
                               Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  for (auto &&elem : blocks) {
    std::visit([&](auto &&block) { block.update_solution(system, y); },
               elem.second);
  }
}

template <typename T>
void Model<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  for (auto &&elem : blocks) {
    std::visit([&](auto &&block) { block.update_constant(system); },
               elem.second);
  }
}

template <typename T>
void Model<T>::update_time(ALGEBRA::SparseSystem<T> &system, T time) {
  for (auto &&elem : blocks) {
    std::visit([&](auto &&block) { block.update_time(system, time); },
               elem.second);
  }
}

template <typename T>
void Model<T>::update_solution(ALGEBRA::SparseSystem<T> &system,
                               Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  for (auto &&elem : blocks) {
    std::visit([&](auto &&block) { block.update_solution(system, y); },
               elem.second);
  }
}

template <typename T>
void Model<T>::to_steady() {
  for (auto &&elem : blocks) {
    std::visit([&](auto &&block) { block.to_steady(); }, elem.second);
  }
}

template <typename T>
void Model<T>::to_unsteady() {
  for (auto &&elem : blocks) {
    std::visit([&](auto &&block) { block.to_unsteady(); }, elem.second);
  }
}

template <typename T>
std::map<std::string, int> Model<T>::get_num_triplets() {
  std::map<std::string, int> num_triplets = {
      {"F", 0},
      {"E", 0},
      {"D", 0},
  };
  for (auto &&elem : blocks) {
    std::visit(
        [&](auto &&block) {
          for (auto &[key, value] : block.num_triplets) {
            num_triplets[key] += value;
          }
        },
        elem.second);
  }
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_MODEL_HPP_