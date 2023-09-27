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
 * @file integrator.hpp
 * @brief ALGEBRA::Integrator source file
 */
#ifndef SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
#define SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_

#include <Eigen/Dense>

#include "../helpers/debug.hpp"
#include "../model/model.hpp"
#include "state.hpp"

namespace ALGEBRA {

/**
 * @brief Generalized-alpha integrator
 *
 * This class handles the time integration scheme for solving 0D blood
 * flow system.
 *
 * Flow rate, pressure, and other hemodynamic quantities in 0D models of
 * vascular anatomies are governed by a system of nonlinear
 * differential-algebraic equations (DAEs):
 *
 * \f[
 * \mathbf{E}(\mathbf{y}, t) \cdot \dot{\mathbf{y}}+\mathbf{F}(\mathbf{y}, t)
 * \cdot \mathbf{y}+\mathbf{c}(\mathbf{y}, t)=\mathbf{0} \f]
 *
 * Here, \f$y\f$ is the vector of solution quantities and \f$\dot{y}\f$ is the
 * time derivative of \f$y\f$. \f$N\f$ is the total number of equations and the
 * total number of global unknowns. The DAE system is solved implicitly using
 * the generalized-\f$\alpha\f$ method \cite JANSEN2000305.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 * `ALGEBRA::SparseSystem`)
 */
template <typename T>
class Integrator {
 private:
  T alpha_m;
  T alpha_f;
  T alpha_m_inv;
  T alpha_f_inv;
  T gamma;
  T gamma_inv;
  T time_step_size;
  T time_step_size_inv;
  T y_dot_coeff;
  T atol;
  T y_init_coeff;
  T ydot_init_coeff;
  int max_iter;
  int size;
  Eigen::Matrix<T, Eigen::Dynamic, 1> y_af;
  Eigen::Matrix<T, Eigen::Dynamic, 1> ydot_am;
  SparseSystem<T> system;
  MODEL::Model<T>* model;

 public:
  /**
   * @brief Construct a new Integrator object
   *
   * @param model The model to simulate
   * @param time_step_size Time step size for generalized-alpha step
   * @param rho Spectral radius for generalized-alpha step
   * @param atol Absolut tolerance for non-linear iteration termination
   * @param max_iter Maximum number of non-linear iterations
   */
  Integrator(MODEL::Model<T>* model, T time_step_size, T rho, T atol,
             int max_iter);

  /**
   * @brief Construct a new Integrator object
   *
   */
  Integrator();

  /**
   * @brief Destroy the Integrator object
   *
   */
  ~Integrator();

  /**
   * @brief Delete dynamically allocated memory (in class member SparseSystem<T>
   * system).
   */
  void clean();

  /**
   * @brief Update integrator parameter and system matrices with model parameter
   * updates.
   *
   * @param time_step_size Time step size for 0D model
   */
  void update_params(T time_step_size);

  /**
   * @brief Perform a time step
   *
   * @param state Current state
   * @param time Current time
   * @return New state
   */
  State<T> step(State<T>& state, T time);
};

template <typename T>
Integrator<T>::Integrator(MODEL::Model<T>* model, T time_step_size, T rho,
                          T atol, int max_iter) {
  this->model = model;
  alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
  alpha_f = 1.0 / (1.0 + rho);
  alpha_m_inv = 1.0 / alpha_m;
  alpha_f_inv = 1.0 / alpha_f;
  gamma = 0.5 + alpha_m - alpha_f;
  gamma_inv = 1.0 / gamma;

  y_init_coeff = alpha_f * 0.5 * time_step_size;
  ydot_init_coeff = (1.0 + alpha_m * ((gamma - 0.5) * gamma_inv - 1.0));

  size = model->dofhandler.size();
  system = SparseSystem<T>(size);
  this->time_step_size = time_step_size;
  this->atol = atol;
  this->max_iter = max_iter;
  time_step_size_inv = 1.0 / time_step_size;

  y_af = Eigen::Matrix<T, Eigen::Dynamic, 1>(size);
  ydot_am = Eigen::Matrix<T, Eigen::Dynamic, 1>(size);

  y_dot_coeff = alpha_m * alpha_f_inv * gamma_inv * time_step_size_inv;

  // Make some memory reservations
  system.reserve(model);
}

template <typename T>
Integrator<T>::Integrator() {}

template <typename T>
Integrator<T>::~Integrator() {}

template <typename T>
void Integrator<T>::clean() {
  // Cannot be in destructor because dynamically allocated pointers will be lost
  // when objects are assigned from temporary objects.
  system.clean();
}

template <typename T>
void Integrator<T>::update_params(T time_step_size) {
  y_init_coeff = alpha_f * 0.5 * time_step_size;
  this->time_step_size = time_step_size;
  time_step_size_inv = 1.0 / time_step_size;
  y_dot_coeff = alpha_m * alpha_f_inv * gamma_inv * time_step_size_inv;
  model->update_constant(system);
  model->update_time(system, 0.0);
}

template <typename T>
State<T> Integrator<T>::step(State<T>& old_state, T time) {
  // Predictor + initiator step
  y_af.setZero();
  ydot_am.setZero();
  y_af += old_state.y + old_state.ydot * y_init_coeff;
  ydot_am += old_state.ydot * ydot_init_coeff;

  // Determine new time
  T new_time = time + alpha_f * time_step_size;

  // Update time-dependent element contributions in system
  model->update_time(system, new_time);

  for (size_t i = 0; i < max_iter; i++) {
    // Update solution-dependent element contribitions
    model->update_solution(system, y_af, ydot_am);

    // Update residuum and check termination criteria
    system.update_residual(y_af, ydot_am);
    if (system.residual.cwiseAbs().maxCoeff() < atol) {
      break;
    }

    // Abort if maximum number of non-linear iterations is reached
    else if (i == max_iter - 1) {
      throw std::runtime_error(
          "Maximum number of non-linear iterations reached.");
    }

    // Determine jacobian
    system.update_jacobian(y_dot_coeff);

    // Solve system
    system.solve();

    // Add increment to solution
    y_af += system.dy;
    ydot_am += system.dy * y_dot_coeff;
  }

  // Set new state
  State<T> new_state = State<T>::Zero(size);
  new_state.y += old_state.y + (y_af - old_state.y) * alpha_f_inv;
  new_state.ydot += old_state.ydot + (ydot_am - old_state.ydot) * alpha_m_inv;

  return new_state;
}
}  // namespace ALGEBRA

#endif  // SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
