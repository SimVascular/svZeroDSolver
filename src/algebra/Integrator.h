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

#include "Model.h"
#include "State.h"

namespace algebra {

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
class Integrator {
 private:
  double alpha_m;
  double alpha_f;
  double alpha_m_inv;
  double alpha_f_inv;
  double gamma;
  double gamma_inv;
  double time_step_size;
  double time_step_size_inv;
  double y_dot_coeff;
  double atol;
  double y_init_coeff;
  double ydot_init_coeff;
  int max_iter;
  int size;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_af;
  Eigen::Matrix<double, Eigen::Dynamic, 1> ydot_am;
  SparseSystem system;
  zd_model::Model* model;

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
  Integrator(zd_model::Model* model, double time_step_size, double rho, double atol,
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
   * @brief Delete dynamically allocated memory (in class member SparseSystem<double>
   * system).
   */
  void clean();

  /**
   * @brief Update integrator parameter and system matrices with model parameter
   * updates.
   *
   * @param time_step_size Time step size for 0D model
   */
  void update_params(double time_step_size);

  /**
   * @brief Perform a time step
   *
   * @param state Current state
   * @param time Current time
   * @return New state
   */
  State step(State& state, double time);
};

}  

#endif  // SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
