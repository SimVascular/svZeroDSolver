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