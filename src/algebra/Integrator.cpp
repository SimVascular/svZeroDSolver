// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "Integrator.h"

Integrator::Integrator(Model* model, double time_step_size, double rho,
                       double atol, int max_iter) {
  this->model = model;
  alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
  alpha_f = 1.0 / (1.0 + rho);
  gamma = 0.5 + alpha_m - alpha_f;
  ydot_init_coeff = 1.0 - 1.0 / gamma;

  y_coeff = gamma * time_step_size;
  y_coeff_jacobian = alpha_f * y_coeff;

  size = model->dofhandler.size();
  system = SparseSystem(size);
  this->time_step_size = time_step_size;
  this->atol = atol;
  this->max_iter = max_iter;

  y_af = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);
  ydot_am = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);

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
  this->time_step_size = time_step_size;
  y_coeff = gamma * time_step_size;
  y_coeff_jacobian = alpha_f * y_coeff;
  model->update_constant(system);
  model->update_time(system, 0.0);
}

State Integrator::step(const State& old_state, double time) {
  // Predictor: Constant y, consistent ydot
  State new_state = State::Zero(size);
  new_state.ydot += old_state.ydot * ydot_init_coeff;
  new_state.y += old_state.y;

  // Determine new time (evaluate terms at generalized mid-point)
  double new_time = time + alpha_f * time_step_size;

  // Evaluate time-dependent element contributions in system
  model->update_time(system, new_time);

  // Count total number of step calls
  n_iter++;

  // Non-linear Newton-Raphson iterations
  for (size_t i = 0; i < max_iter; i++) {
    // Initiator: Evaluate the iterates at the intermediate time levels
    ydot_am.setZero();
    y_af.setZero();
    ydot_am += old_state.ydot + (new_state.ydot - old_state.ydot) * alpha_m;
    y_af += old_state.y + (new_state.y - old_state.y) * alpha_f;

    // Update solution-dependent element contribitions
    model->update_solution(system, y_af, ydot_am);

    // Evaluate residual
    system.update_residual(y_af, ydot_am);

    // Check termination criterium
    if (system.residual.cwiseAbs().maxCoeff() < atol) {
      break;
    }

    // Abort if maximum number of non-linear iterations is reached
    else if (i == max_iter - 1) {
      throw std::runtime_error(
          "Maximum number of non-linear iterations reached.");
    }

    // Evaluate Jacobian
    system.update_jacobian(alpha_m, y_coeff_jacobian);

    // Solve system for increment in ydot
    system.solve();

    // Perform post-solve actions on blocks
    model->post_solve(new_state.y);

    // Update the solution
    new_state.ydot += system.dydot;
    new_state.y += system.dydot * y_coeff;

    // Count total number of nonlinear iterations
    n_nonlin_iter++;
  }

  return new_state;
}

double Integrator::avg_nonlin_iter() {
  return (double)n_nonlin_iter / (double)n_iter;
}
