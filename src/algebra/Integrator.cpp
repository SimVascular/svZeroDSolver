// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "Integrator.h"

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

  size = model->dofhandler.size();
  system = SparseSystem(size);
  this->time_step_size = time_step_size;
  this->model->set_time_step_size(time_step_size);
  this->atol = atol;
  this->max_iter = max_iter;
  this->max_iter_error_to_warning = max_iter_error_to_warning;

  y_af = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);
  ydot_am = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);

  // The convergence check skips equations that only define a diagnostic
  // output, see DOFHandler::register_equation.
  const auto& is_output = model->dofhandler.equation_is_output;
  residual_check_dofs.reserve(size);
  for (int i = 0; i < size; i++) {
    if (i >= (int)is_output.size() || !is_output[i]) {
      residual_check_dofs.push_back(i);
    }
  }

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

double Integrator::evaluate_residual(const State& st, const State& old_state,
                                     int& res_dof) {
  // Initiator: Evaluate the iterates at the intermediate time levels
  ydot_am.setZero();
  y_af.setZero();
  ydot_am += old_state.ydot + (st.ydot - old_state.ydot) * alpha_m;
  y_af += old_state.y + (st.y - old_state.y) * alpha_f;

  // Update solution-dependent element contributions
  model->update_solution(system, y_af, ydot_am);

  // Evaluate residual
  system.update_residual(y_af, ydot_am);

  double res = 0.0;
  res_dof = residual_check_dofs.empty() ? 0 : residual_check_dofs.front();
  for (int i : residual_check_dofs) {
    double abs_res = std::abs(system.residual[i]);
    if (abs_res > res) {
      res = abs_res;
      res_dof = i;
    }
  }
  return res;
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

  // Per-step convergence history
  int iter_this_step = 0;
  double res_this_step = 0.0;
  int res_dof_this_step = 0;
  std::vector<double> history_this_step;

  // The residual is normally evaluated at the top of each iteration. The line
  // search already evaluates it at the accepted iterate, so it can be reused.
  bool have_residual = false;

  // Non-linear Newton-Raphson iterations
  for (int i = 0; i < max_iter; i++) {
    if (!have_residual) {
      res_this_step =
          evaluate_residual(new_state, old_state, res_dof_this_step);
    }
    have_residual = false;

    if (record_history) {
      history_this_step.push_back(res_this_step);
    }

    // Check termination criterium
    if (res_this_step < atol) {
      break;
    }

    // Abort if maximum number of non-linear iterations is reached
    else if (i == max_iter - 1) {
      if (max_iter_error_to_warning) {
        std::cout << "Warning: Maximum number of non-linear iterations reached "
                     "at time "
                  << time << std::endl;
      } else {
        throw std::runtime_error(
            "Maximum number of non-linear iterations reached at time " +
            std::to_string(time));
      }
    }

    // Evaluate Jacobian
    system.update_jacobian(alpha_m, y_coeff_jacobian);

    // Solve system for increment in ydot
    system.solve();

    // Perform post-solve actions on blocks
    model->post_solve(new_state.y);

    if (max_line_search_steps > 0) {
      // Affine-invariant (natural) line search after Deuflhard: shorten the
      // Newton step until the simplified Newton correction
      // J(y_k)^-1 r(y_k + lambda dy) is shorter than the full step. Measuring
      // progress in the correction rather than in the residual makes the test
      // independent of how the individual equations are scaled, which matters
      // here because pressure and flow residuals differ by many orders of
      // magnitude. The full step is tried first, so a well-behaved iteration
      // only pays for one triangular solve.
      const double step_norm = system.dydot.norm();
      Eigen::Matrix<double, Eigen::Dynamic, 1> full_step = system.dydot;
      double lambda = 1.0;
      State trial = new_state;
      for (int k = 0; k <= max_line_search_steps; k++) {
        trial.ydot = new_state.ydot + full_step * lambda;
        trial.y = new_state.y + full_step * (y_coeff * lambda);
        int trial_dof = 0;
        double trial_res = evaluate_residual(trial, old_state, trial_dof);
        double corr_norm =
            system.apply_inverse_jacobian(system.residual).norm();
        if (corr_norm <= line_search_factor * step_norm ||
            k == max_line_search_steps || step_norm == 0.0) {
          new_state = trial;
          res_this_step = trial_res;
          res_dof_this_step = trial_dof;
          have_residual = true;
          break;
        }
        lambda *= 0.5;
      }
    } else {
      // Update the solution
      new_state.ydot += system.dydot;
      new_state.y += system.dydot * y_coeff;
    }

    // Count total number of nonlinear iterations
    n_nonlin_iter++;
    iter_this_step++;
  }

  // Record the per-step convergence history
  time_per_step.push_back(time);
  iter_per_step.push_back(iter_this_step);
  res_per_step.push_back(res_this_step);
  res_dof_per_step.push_back(res_dof_this_step);
  if (record_history) {
    residual_history.push_back(std::move(history_this_step));
  }

  model->accept_timestep(new_state.y);
  return new_state;
}

double Integrator::avg_nonlin_iter() {
  return (double)n_nonlin_iter / (double)n_iter;
}
