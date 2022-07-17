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
 * @file densesystem.hpp
 * @brief ALGEBRA::DenseSystem source file
 */
#ifndef SVZERODSOLVER_ALGREBRA_DENSESYSTEM_HPP_
#define SVZERODSOLVER_ALGREBRA_DENSESYSTEM_HPP_

#include <Eigen/Dense>

// Forward declaration of Model
namespace MODEL {
template <typename T>
class Model;
}

namespace ALGEBRA {

/**
 * @brief Dense system
 *
 * This class contains all attributes and methods to create, modify, and
 * solve dense systems.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class DenseSystem {
 public:
  /**
   * @brief Construct a new Dense System object
   *
   */
  DenseSystem();

  /**
   * @brief Construct a new Dense System object
   *
   * @param n Size of the system
   */
  DenseSystem(unsigned int n);

  /**
   * @brief Destroy the Dense System object
   *
   */
  ~DenseSystem();

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> F;  ///< System matrix F
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> E;  ///< System matrix E
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> D;  ///< System matrix D
  Eigen::Matrix<T, Eigen::Dynamic, 1> C;               ///< System matrix C

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
      jacobian;                                  ///< Jacobian of the system
  Eigen::Matrix<T, Eigen::Dynamic, 1> residual;  ///< Residual of the system
  Eigen::Matrix<T, Eigen::Dynamic, 1> dy;  ///< Solution increment of the system

  /**
   * @brief Reserve memory in system matrices based on number of triplets
   *
   * @param model The model to reserve space for in the system
   */
  void reserve(MODEL::Model<T> &model);

  /**
   * @brief Update the residual of the system
   *
   * @param y Vector of current solution quantities
   * @param ydot Derivate of y
   */
  void update_residual(Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &ydot);

  /**
   * @brief Update the jacobian of the system
   *
   * @param e_coeff Coefficent for system matrix \ref E
   */
  void update_jacobian(T e_coeff);

  /**
   * @brief Solve the system
   */
  void solve();
};

template <typename T>
DenseSystem<T>::DenseSystem() {}

template <typename T>
DenseSystem<T>::DenseSystem(unsigned int n) {
  F = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
  E = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
  D = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
  C = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);

  jacobian = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
  residual = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
  dy = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
}

template <typename T>
DenseSystem<T>::~DenseSystem() {}

template <typename T>
void DenseSystem<T>::reserve(MODEL::Model<T> &model) {
  model.update_constant(*this);
  model.update_time(*this, 0.0);
  Eigen::Matrix<T, Eigen::Dynamic, 1> dummy_y =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Ones(residual.size());
  model.update_solution(*this, dummy_y);
}

template <typename T>
void DenseSystem<T>::update_residual(
    Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &ydot) {
  residual = -(E * ydot) - (F * y) - C;
}

template <typename T>
void DenseSystem<T>::update_jacobian(T e_coeff) {
  jacobian = F + D + E * e_coeff;
}

template <typename T>
void DenseSystem<T>::solve() {
  dy = jacobian.colPivHouseholderQr().solve(residual);
}

}  // namespace ALGEBRA

#include "../model/model.hpp"

#endif  // SVZERODSOLVER_ALGREBRA_DENSESYSTEM_HPP_
