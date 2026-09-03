// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ImpedanceBC.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "Model.h"

void ImpedanceBC::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 1, {});
}

void ImpedanceBC::update_constant(SparseSystem& system,
                                  std::vector<double>& parameters) {
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
}

void ImpedanceBC::update_time(SparseSystem& system,
                              std::vector<double>& parameters) {
  ensure_initialized();

  double conv_sum = 0.0;
  // Match Olufsen startup semantics: do not evaluate lagged terms until one
  // full cycle of accepted data has elapsed.
  if (num_accepted_steps >= num_period_steps) {
    for (int m = 1; m < num_kernel_terms; m++) {
      conv_sum += kernel[m] * lagged_flow(m);
    }
  }

  // Kernel entries are discrete-time convolution coefficients (Olufsen
  // periodic kernel), so no extra dt scaling is applied in the algebraic BC.
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -kernel[0];
  system.C(global_eqn_ids[0]) = -(pd + conv_sum);
}

void ImpedanceBC::accept_timestep(
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y) {
  ensure_initialized();

  const double flow_accepted = y[global_var_ids[1]];
  head = (head + 1) % num_period_steps;
  flow_history[head] = flow_accepted;
  committed_samples = std::min(committed_samples + 1, num_period_steps);
  num_accepted_steps += 1;
}

std::vector<double> ImpedanceBC::get_persistent_state() const {
  if (!initialized) {
    return {};
  }

  std::vector<double> state;
  state.reserve(6 + flow_history.size());
  state.push_back(dt);
  state.push_back(static_cast<double>(num_period_steps));
  state.push_back(static_cast<double>(num_kernel_terms));
  state.push_back(static_cast<double>(head));
  state.push_back(static_cast<double>(committed_samples));
  state.push_back(static_cast<double>(num_accepted_steps));
  state.insert(state.end(), flow_history.begin(), flow_history.end());
  return state;
}

void ImpedanceBC::set_persistent_state(const std::vector<double>& state) {
  if (state.empty()) {
    initialized = false;
    dt = -1.0;
    num_period_steps = -1;
    num_kernel_terms = -1;
    head = -1;
    committed_samples = 0;
    num_accepted_steps = 0;
    flow_history.clear();
    return;
  }

  if (state.size() < 5) {
    throw std::runtime_error("Invalid persistent state for IMPEDANCE block " +
                             get_name() + ".");
  }

  const double new_dt = state[0];
  const int new_num_period_steps = static_cast<int>(std::llround(state[1]));
  const int new_num_kernel_terms = static_cast<int>(std::llround(state[2]));
  const int new_head = static_cast<int>(std::llround(state[3]));
  const int new_committed_samples = static_cast<int>(std::llround(state[4]));
  const bool has_num_accepted_steps =
      (state.size() == static_cast<size_t>(6 + new_num_period_steps));
  const long long new_num_accepted_steps =
      has_num_accepted_steps
          ? static_cast<long long>(std::llround(state[5]))
          : static_cast<long long>(new_committed_samples);

  if ((new_dt <= 0.0) || (new_num_period_steps <= 0) ||
      (new_num_kernel_terms <= 0) ||
      (new_num_kernel_terms > new_num_period_steps) || (new_head < 0) ||
      (new_head >= new_num_period_steps) || (new_committed_samples < 0) ||
      (new_committed_samples > new_num_period_steps) ||
      (new_num_accepted_steps < 0)) {
    throw std::runtime_error("Corrupted persistent state for IMPEDANCE block " +
                             get_name() + ".");
  }

  const size_t expected_size = has_num_accepted_steps
                                   ? static_cast<size_t>(6 + new_num_period_steps)
                                   : static_cast<size_t>(5 + new_num_period_steps);
  if (state.size() != expected_size) {
    throw std::runtime_error("Persistent state size mismatch for IMPEDANCE "
                             "block " +
                             get_name() + ".");
  }

  if (configured && (kernel.size() != static_cast<size_t>(new_num_period_steps))) {
    throw std::runtime_error(
        "Persistent IMPEDANCE state is incompatible with configured kernel "
        "size in block " +
        get_name() + ".");
  }

  dt = new_dt;
  num_period_steps = new_num_period_steps;
  num_kernel_terms = new_num_kernel_terms;
  head = new_head;
  committed_samples = new_committed_samples;
  num_accepted_steps = new_num_accepted_steps;
  flow_history.assign(state.begin() + (has_num_accepted_steps ? 6 : 5),
                      state.end());
  initialized = true;
}

void ImpedanceBC::configure(const std::vector<double>& z, double pd,
                            const std::string& convolution_mode,
                            int num_kernel_terms) {
  if (z.empty()) {
    throw std::runtime_error("IMPEDANCE block " + get_name() +
                             " requires non-empty kernel `z`.");
  }
  for (size_t i = 0; i < z.size(); i++) {
    if (!std::isfinite(z[i])) {
      throw std::runtime_error("IMPEDANCE block " + get_name() +
                               " has non-finite kernel value at index " +
                               std::to_string(i) + ".");
    }
  }

  std::string mode_lower = convolution_mode;
  std::transform(mode_lower.begin(), mode_lower.end(), mode_lower.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (mode_lower.empty()) {
    mode_lower = "exact";
  }

  if (mode_lower == "exact") {
    mode = ConvolutionMode::exact;
  } else if (mode_lower == "truncated") {
    mode = ConvolutionMode::truncated;
  } else {
    throw std::runtime_error("IMPEDANCE block " + get_name() +
                             " has invalid `convolution_mode` value `" +
                             convolution_mode + "`. Supported values are "
                                                "`exact` and `truncated`.");
  }

  if ((mode == ConvolutionMode::truncated) && (num_kernel_terms <= 0)) {
    throw std::runtime_error("IMPEDANCE block " + get_name() +
                             " in truncated mode requires "
                             "`num_kernel_terms > 0`.");
  }

  kernel = z;
  this->pd = pd;
  num_kernel_terms_input = num_kernel_terms;
  configured = true;

  // Reset runtime state; it will be re-initialized using current dt.
  initialized = false;
  dt = -1.0;
  num_period_steps = -1;
  this->num_kernel_terms = -1;
  head = -1;
  committed_samples = 0;
  num_accepted_steps = 0;
  flow_history.clear();
}

void ImpedanceBC::ensure_initialized() {
  if (!configured) {
    throw std::runtime_error("IMPEDANCE block " + get_name() +
                             " was not configured. Missing or invalid "
                             "impedance kernel input.");
  }

  const double model_dt = model->get_time_step_size();
  if (model_dt <= 0.0) {
    throw std::runtime_error("IMPEDANCE block " + get_name() +
                             " requires a positive solver time-step size.");
  }
  const double period = model->cardiac_cycle_period;
  if (period <= 0.0) {
    throw std::runtime_error("IMPEDANCE block " + get_name() +
                             " requires `simulation_parameters.cardiac_period "
                             "> 0`.");
  }

  if (!initialized) {
    const double period_over_dt = period / model_dt;
    const long rounded = std::lround(period_over_dt);
    if (rounded <= 0) {
      throw std::runtime_error("IMPEDANCE block " + get_name() +
                               " has invalid period-to-dt ratio.");
    }

    const double mismatch =
        std::abs(period_over_dt - static_cast<double>(rounded));
    if (mismatch > 1.0e-8) {
      throw std::runtime_error("IMPEDANCE block " + get_name() +
                               " requires period/dt to be an integer. Got "
                               "period=" +
                               std::to_string(period) +
                               ", dt=" + std::to_string(model_dt) + ".");
    }

    const int dt_period_steps = static_cast<int>(rounded);
    if (model->impedance_points_per_cycle > 0) {
      const int configured_period_steps = model->impedance_points_per_cycle - 1;
      if (configured_period_steps <= 0) {
        throw std::runtime_error(
            "IMPEDANCE block " + get_name() +
            " requires `simulation_parameters.number_of_time_pts_per_cardiac_cycle > 1`.");
      }
      if (configured_period_steps != dt_period_steps) {
        throw std::runtime_error(
            "IMPEDANCE block " + get_name() +
            " requires `simulation_parameters.number_of_time_pts_per_cardiac_cycle - 1` "
            "to match cardiac_period / dt. Expected " +
            std::to_string(dt_period_steps + 1) + ", got " +
            std::to_string(model->impedance_points_per_cycle) + ".");
      }
      num_period_steps = configured_period_steps;
    } else {
      num_period_steps = dt_period_steps;
    }
    if (kernel.size() != static_cast<size_t>(num_period_steps)) {
      throw std::runtime_error(
          "IMPEDANCE block " + get_name() +
          " requires kernel length `z.size()` to match one period in time "
          "steps. Expected " +
          std::to_string(num_period_steps) + ", got " +
          std::to_string(kernel.size()) + ".");
    }

    dt = model_dt;
    if (mode == ConvolutionMode::exact) {
      num_kernel_terms = num_period_steps;
    } else {
      num_kernel_terms = num_kernel_terms_input;
    }
    if ((num_kernel_terms <= 0) || (num_kernel_terms > num_period_steps)) {
      throw std::runtime_error(
          "IMPEDANCE block " + get_name() +
          " has invalid `num_kernel_terms`. It must be in [1, " +
          std::to_string(num_period_steps) + "].");
    }

    flow_history.assign(num_period_steps, 0.0);
    head = num_period_steps - 1;
    committed_samples = 0;
    num_accepted_steps = 0;
    initialized = true;
    return;
  }

  const double denom =
      std::max({1.0, std::abs(model_dt), std::abs(dt)});
  const double rel_change = std::abs(model_dt - dt) / denom;
  if (rel_change > 1.0e-12) {
    throw std::runtime_error(
        "IMPEDANCE block " + get_name() +
        " does not support changing dt after initialization. "
        "Reinitialize the model when external step size changes.");
  }
}

double ImpedanceBC::lagged_flow(int lag) const {
  if (lag <= 0) {
    throw std::runtime_error("Internal error: lag must be positive in "
                             "IMPEDANCE block " +
                             get_name() + ".");
  }
  if (lag > committed_samples) {
    return 0.0;
  }

  int idx = (head - (lag - 1)) % num_period_steps;
  if (idx < 0) {
    idx += num_period_steps;
  }
  return flow_history[idx];
}
