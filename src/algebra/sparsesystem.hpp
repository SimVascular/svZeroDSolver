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
 * @file sparsesystem.hpp
 * @brief ALGEBRA::SparseSystem source file
 */
#ifndef SVZERODSOLVER_ALGREBRA_SPARSESYSTEM_HPP_
#define SVZERODSOLVER_ALGREBRA_SPARSESYSTEM_HPP_

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>

// Forward declaration of Model
namespace MODEL {
template <typename T>
class Model;
}

namespace ALGEBRA {
/**
 * @brief Sparse system
 *
 * This class contains all attributes and methods to create, modify, and
 * solve sparse systems.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class SparseSystem {
 public:
  /**
   * @brief Construct a new Sparse System object
   *
   */
  SparseSystem();

  /**
   * @brief Construct a new Sparse System object
   *
   * @param n Size of the system
   */
  SparseSystem(unsigned int n);

  /**
   * @brief Destroy the Sparse System object
   *
   */
  ~SparseSystem();

  Eigen::SparseMatrix<T> F;               ///< System matrix F
  Eigen::SparseMatrix<T> E;               ///< System matrix E
  Eigen::SparseMatrix<T> D;               ///< System matrix D
  Eigen::Matrix<T, Eigen::Dynamic, 1> C;  ///< System vector C

  Eigen::SparseMatrix<T> jacobian;               ///< Jacobian of the system
  Eigen::Matrix<T, Eigen::Dynamic, 1> residual;  ///< Residual of the system
  Eigen::Matrix<T, Eigen::Dynamic, 1> dy;  ///< Solution increment of the system

  Eigen::SparseLU<Eigen::SparseMatrix<T>> *solver =
      new Eigen::SparseLU<Eigen::SparseMatrix<T>>();  ///< Linear solver

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

  /**
   * @brief Delete dynamically allocated memory (class member
   * Eigen::SparseLU<Eigen::SparseMatrix<T>> *solver)
   */
  void clean();
};

template <typename T>
SparseSystem<T>::SparseSystem() {}

template <typename T>
SparseSystem<T>::SparseSystem(unsigned int n) {
  F = Eigen::SparseMatrix<T>(n, n);
  E = Eigen::SparseMatrix<T>(n, n);
  D = Eigen::SparseMatrix<T>(n, n);
  C = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);

  jacobian = Eigen::SparseMatrix<T>(n, n);
  residual = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
  dy = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
}

template <typename T>
SparseSystem<T>::~SparseSystem() {}

template <typename T>
void SparseSystem<T>::clean() {
  // Cannot be in destructor because dynamically allocated pointers will be lost
  // when objects are assigned from temporary objects.
  delete solver;
}

template <typename T>
void SparseSystem<T>::reserve(MODEL::Model<T> &model) {
  auto num_triplets = model.get_num_triplets();
  F.reserve(num_triplets["F"]);
  E.reserve(num_triplets["E"]);
  D.reserve(num_triplets["D"]);
  model.update_constant(*this);
  model.update_time(*this, 0.0);
  Eigen::Matrix<T, Eigen::Dynamic, 1> dummy_y =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Ones(residual.size());
  model.update_solution(*this, dummy_y);
  F.makeCompressed();
  E.makeCompressed();
  D.makeCompressed();
  jacobian.reserve(num_triplets["F"] + num_triplets["E"]);  // Just an estimate
  update_jacobian(1.0);  // Update it once to have sparsity pattern
  jacobian.makeCompressed();
  solver->analyzePattern(jacobian);  // Let solver analyze pattern
}

template <typename T>
void SparseSystem<T>::update_residual(
    Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &ydot) {
  residual.setZero();
  residual -= C;
  residual.noalias() -= E * ydot;
  residual.noalias() -= F * y;
}

template <typename T>
void SparseSystem<T>::update_jacobian(T e_coeff) {
  jacobian.setZero();
  jacobian += F + D + E * e_coeff;
}

template <typename T>
void SparseSystem<T>::solve() {
  solver->factorize(jacobian);
  dy.setZero();
  dy += solver->solve(residual);
}

}  // namespace ALGEBRA

#include "../model/model.hpp"

#endif  // SVZERODSOLVER_ALGREBRA_SPARSESYSTEM_HPP_
