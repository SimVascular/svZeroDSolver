// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "Integrator.h"

// ---------------------------------------------------------------------------
// TimeIntegrator (shared machinery)
// ---------------------------------------------------------------------------

void TimeIntegrator::init(Model* model, double time_step_size, double rho,
                          double atol, int max_iter,
                          bool max_iter_error_to_warning) {
  this->model = model;
  alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
  alpha_f = 1.0 / (1.0 + rho);
  gamma = 0.5 + alpha_m - alpha_f;
  ydot_init_coeff = 1.0 - 1.0 / gamma;

  size = model->dofhandler.size();
  system = SparseSystem(size);
  this->time_step_size = time_step_size;
  this->model->set_time_step_size(time_step_size);
  this->atol = atol;
  this->max_iter = max_iter;
  this->max_iter_error_to_warning = max_iter_error_to_warning;

  y_af = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);
  ydot_am = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);

  system.reserve(model);
}

void TimeIntegrator::clean() { system.clean(); }

void TimeIntegrator::update_params(double time_step_size) {
  this->time_step_size = time_step_size;
  model->set_time_step_size(time_step_size);
  model->update_constant(system);
  model->update_time(system, 0.0);
}

double TimeIntegrator::avg_nonlin_iter() {
  return (double)n_nonlin_iter / (double)n_iter;
}

double TimeIntegrator::assemble_residual(const State& old_state,
                                         const State& new_state) {
  ydot_am.setZero();
  y_af.setZero();
  ydot_am += old_state.ydot + (new_state.ydot - old_state.ydot) * alpha_m;
  y_af += old_state.y + (new_state.y - old_state.y) * alpha_f;
  model->update_solution(system, y_af, ydot_am);
  system.update_residual(y_af, ydot_am);
  return system.residual.cwiseAbs().maxCoeff();
}

// ---------------------------------------------------------------------------
// GeneralizedAlpha: the generalized-alpha scheme (\cite JANSEN2000305), the
// default svZeroDSolver integrator.
// ---------------------------------------------------------------------------

GeneralizedAlpha::GeneralizedAlpha(Model* model, double time_step_size,
                                   double rho, double atol, int max_iter,
                                   bool max_iter_error_to_warning) {
  init(model, time_step_size, rho, atol, max_iter, max_iter_error_to_warning);
}

State GeneralizedAlpha::step(const State& old_state, double time) {
  // Step-size-dependent coefficients (map the ydot increment to y)
  const double y_coeff = gamma * time_step_size;
  const double y_coeff_jacobian = alpha_f * y_coeff;

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
    // Evaluate residual at the intermediate time levels
    double res_norm = assemble_residual(old_state, new_state);

    // Check termination criterium
    if (res_norm < atol) {
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

    // Update the solution
    new_state.ydot += system.dydot;
    new_state.y += system.dydot * y_coeff;

    // Count total number of nonlinear iterations
    n_nonlin_iter++;
  }
  model->accept_timestep(new_state.y);
  return new_state;
}

// ---------------------------------------------------------------------------
// ConsistentStiffIntegrator (damped Newton + adaptive sub-stepping)
// ---------------------------------------------------------------------------

ConsistentStiffIntegrator::ConsistentStiffIntegrator(
    Model* model, double time_step_size, double rho, double atol, int max_iter,
    bool max_iter_error_to_warning) {
  // Spectral radius 0 (L-stable) regardless of rho: the stiff active law of
  // \cite genet23 needs the high-frequency damping.
  (void)rho;
  init(model, time_step_size, 0.0, atol, max_iter, max_iter_error_to_warning);
}

bool ConsistentStiffIntegrator::damped_step(const State& old_state,
                                            State& new_state, double dt,
                                            double time) {
  // Step-size-dependent coefficients of this (sub-)step
  const double y_coeff = gamma * dt;
  const double y_coeff_jacobian = alpha_f * y_coeff;

  // Predictor: constant y, consistent ydot.
  new_state = State::Zero(size);
  new_state.ydot += old_state.ydot * ydot_init_coeff;
  new_state.y += old_state.y;

  model->update_time(system, time + alpha_f * dt);

  for (int i = 0; i < max_iter; i++) {
    double res_norm = assemble_residual(old_state, new_state);
    if (res_norm < atol) {
      return true;
    }

    system.update_jacobian(alpha_m, y_coeff_jacobian);
    system.solve();
    model->post_solve(new_state.y);
    n_nonlin_iter++;

    // Backtracking line search on the residual infinity-norm (damped Newton).
    double lambda = 1.0;
    bool accepted = false;
    for (int bt = 0; bt < max_backtrack; bt++) {
      State trial = new_state;
      trial.ydot += system.dydot * lambda;
      trial.y += system.dydot * (y_coeff * lambda);
      double trial_norm = assemble_residual(old_state, trial);
      if (trial_norm < res_norm) {
        new_state = trial;
        accepted = true;
        break;
      }
      lambda *= 0.5;
    }
    // If no damping reduced the residual, take the undamped step and continue.
    if (!accepted) {
      new_state.ydot += system.dydot;
      new_state.y += system.dydot * y_coeff;
    }
  }
  return false;
}

State ConsistentStiffIntegrator::step(const State& old_state, double time) {
  n_iter++;
  int nsub = 1;
  for (int level = 0; level <= max_subdiv; level++, nsub *= 2) {
    const double dt = time_step_size / nsub;
    State st = old_state;
    bool ok = true;
    for (int k = 0; k < nsub; k++) {
      State next;
      if (!damped_step(st, next, dt, time + k * dt)) {
        ok = false;
        break;
      }
      st = next;
    }
    if (ok) {
      // Blocks with persistent memory (e.g. ImpedanceBC) see one accepted step
      // of the nominal size; they assume no sub-stepping.
      model->accept_timestep(st.y);
      return st;
    }
  }
  throw std::runtime_error(
      "Consistent-stiff integrator failed to converge (max sub-steps) at "
      "time " +
      std::to_string(time));
}

// ---------------------------------------------------------------------------
// Factory
// ---------------------------------------------------------------------------

std::unique_ptr<TimeIntegrator> make_integrator(
    IntegratorType type, Model* model, double time_step_size, double rho,
    double atol, int max_iter, bool max_iter_error_to_warning) {
  switch (type) {
    case IntegratorType::consistent_stiff:
      return std::make_unique<ConsistentStiffIntegrator>(
          model, time_step_size, rho, atol, max_iter,
          max_iter_error_to_warning);
    case IntegratorType::generalized_alpha:
    default:
      return std::make_unique<GeneralizedAlpha>(model, time_step_size, rho,
                                                atol, max_iter,
                                                max_iter_error_to_warning);
  }
}
