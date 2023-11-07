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

#include "Integrator.h"

Integrator::Integrator(Model* model, double time_step_size, double rho,
                       double atol, int max_iter) {
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
  system = SparseSystem(size);
  this->time_step_size = time_step_size;
  this->atol = atol;
  this->max_iter = max_iter;
  time_step_size_inv = 1.0 / time_step_size;

  y_af = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);
  ydot_am = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);

  y_dot_coeff = alpha_m * alpha_f_inv * gamma_inv * time_step_size_inv;

  // Make some memory reservations
  system.reserve(model);
}

// Must declare default constructord and dedtructor
// because of Eigen.
Integrator::Integrator() {}
Integrator::~Integrator() {}

void Integrator::clean() {
  // Cannot be in destructor because dynamically allocated pointers will be lost
  // when objects are assigned from temporary objects.
  system.clean();
}

void Integrator::update_params(double time_step_size) {
  y_init_coeff = alpha_f * 0.5 * time_step_size;
  this->time_step_size = time_step_size;
  time_step_size_inv = 1.0 / time_step_size;
  y_dot_coeff = alpha_m * alpha_f_inv * gamma_inv * time_step_size_inv;
  model->update_constant(system);
  model->update_time(system, 0.0);
}

State Integrator::step(const State& old_state, double time) {
  // Predictor + initiator step
  y_af.setZero();
  ydot_am.setZero();
  y_af += old_state.y + old_state.ydot * y_init_coeff;
  ydot_am += old_state.ydot * ydot_init_coeff;

  // Determine new time
  double new_time = time + alpha_f * time_step_size;

  // Update time-dependent element contributions in system
  model->update_time(system, new_time);

  // Count total number of step calls
  n_iter++;

  for (size_t i = 0; i < max_iter; i++) {
    // Count total number of nonlinear iterations
    n_nonlin_iter++;

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
  State new_state = State::Zero(size);
  new_state.y += old_state.y + (y_af - old_state.y) * alpha_f_inv;
  new_state.ydot += old_state.ydot + (ydot_am - old_state.ydot) * alpha_m_inv;

  return new_state;
}

double Integrator::avg_nonlin_iter() {
  return (double)n_nonlin_iter / (double)n_iter;
}
