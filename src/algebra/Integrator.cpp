// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "Integrator.h"

#include <sstream>

Integrator::Integrator(Model* model, double time_step_size, double rho,
                       double atol, int max_iter,
                       bool max_iter_error_to_warning) {
  this->model = model;
  alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
  alpha_f = 1.0 / (1.0 + rho);
  gamma = 0.5 + alpha_m - alpha_f;
  ydot_init_coeff = 1.0 - 1.0 / gamma;

  y_coeff = gamma * time_step_size;
  y_coeff_jacobian = alpha_f * y_coeff;

  // The system holds only the solved variables, in the leading block of y
  size = model->dofhandler.get_num_variables();
  num_equations = model->dofhandler.get_num_equations();
  system = SparseSystem(num_equations);
  this->time_step_size = time_step_size;
  this->model->set_time_step_size(time_step_size);
  this->atol = atol;
  this->max_iter = max_iter;
  this->max_iter_error_to_warning = max_iter_error_to_warning;

  y_af = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);
  ydot_am = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);

  // Make some memory reservations
  system.reserve(model);
}

// Must declare default constructors and dedtructor
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
  model->set_time_step_size(time_step_size);
  y_coeff = gamma * time_step_size;
  y_coeff_jacobian = alpha_f * y_coeff;
  model->update_constant(system);
  model->update_time(system, 0.0);
}

void Integrator::evaluate_residual(const State& old_state,
                                   const State& new_state) {
  // Initiator: Evaluate the iterates at the intermediate time levels
  ydot_am.setZero();
  y_af.setZero();
  ydot_am += old_state.ydot + (new_state.ydot - old_state.ydot) * alpha_m;
  y_af += old_state.y + (new_state.y - old_state.y) * alpha_f;

  // Update solution-dependent element contribitions
  model->update_solution(system, y_af, ydot_am);

  // Evaluate residual
  system.update_residual(y_af, ydot_am);
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

  // Convergence history of this step
  int iter_this_step = 0;
  double res_this_step = 0.0;
  Eigen::Index res_dof_this_step = 0;
  bool converged = false;

  // Non-linear Newton-Raphson iterations
  for (int i = 0; i < max_iter; i++) {
    // Evaluate the residual of the current iterate
    evaluate_residual(old_state, new_state);
    res_this_step = system.residual.cwiseAbs().maxCoeff(&res_dof_this_step);

    // Check termination criterium
    if (res_this_step < atol) {
      converged = true;
      break;
    }

    // Evaluate Jacobian
    system.update_jacobian(alpha_m, y_coeff_jacobian);

    // Solve system for increment in ydot
    system.solve();

    // Perform post-solve actions on blocks
    model->post_solve(new_state.y);

    // Update the solution
    new_state.ydot.head(num_equations) += system.dydot;
    new_state.y.head(num_equations) += system.dydot * y_coeff;

    // Count total number of nonlinear iterations
    n_nonlin_iter++;
    iter_this_step++;
  }

  // The loop above leaves the update of its last iteration unchecked. Evaluate
  // the residual once more so that the recorded residual always belongs to the
  // state that is returned, and so that a step that converged with that last
  // update is not reported as a failure.
  if (!converged) {
    evaluate_residual(old_state, new_state);
    res_this_step = system.residual.cwiseAbs().maxCoeff(&res_dof_this_step);
    converged = res_this_step < atol;
  }

  // Record how hard this step was to solve
  history.push_back(time, iter_this_step, res_this_step,
                    static_cast<int>(res_dof_this_step));

  // Abort if the maximum number of non-linear iterations was not enough
  if (!converged) {
    std::ostringstream message;
    message << "Maximum number of non-linear iterations reached at time "
            << time << " (residual " << std::scientific << res_this_step << ")";

    if (max_iter_error_to_warning) {
      std::cout << "Warning: " << message.str() << std::endl;
    } else {
      throw std::runtime_error(message.str());
    }
  }

  // Fill the post-processing variables from the converged state
  model->update_outputs(new_state.y, new_state.ydot);

  model->accept_timestep(new_state.y);
  return new_state;
}

double Integrator::avg_nonlin_iter() {
  return (double)n_nonlin_iter / (double)n_iter;
}
