#include "Solver.h"

#include "csv_writer.h"

Solver::Solver(const nlohmann::json& config) {
  DEBUG_MSG("Read simulation parameters");
  simparams = load_simulation_params(config);
  DEBUG_MSG("Load model");
  model = Model();
  load_simulation_model(config, model);
  DEBUG_MSG("Load initial condition");
  initial_state = load_initial_condition(config, model);

  DEBUG_MSG("Cardiac cycle period " << model.cardiac_cycle_period);

  // Calculate time step size
  if (!simparams.sim_coupled) {
    simparams.sim_time_step_size = model.cardiac_cycle_period /
                                   (double(simparams.sim_pts_per_cycle) - 1.0);
  } else {
    simparams.sim_time_step_size = simparams.sim_external_step_size /
                                   (double(simparams.sim_num_time_steps) - 1.0);
  }
  sanity_checks();
}

Solver::~Solver() {}

void Solver::run() {
  auto state = initial_state;

  // Create steady initial
  if (simparams.sim_steady_initial) {
    DEBUG_MSG("Calculate steady initial condition");
    double time_step_size_steady = model.cardiac_cycle_period / 10.0;
    model.to_steady();

    Integrator integrator_steady(&model, time_step_size_steady, 0.1,
                                 simparams.sim_abs_tol, simparams.sim_nliter);

    for (int i = 0; i < 31; i++) {
      state = integrator_steady.step(state, time_step_size_steady * double(i));
    }

    model.to_unsteady();
  }

  // Set-up integrator
  DEBUG_MSG("Setup time integration");
  Integrator integrator(&model, simparams.sim_time_step_size, 0.1,
                        simparams.sim_abs_tol, simparams.sim_nliter);

  // Initialize loop
  states = std::vector<State>();
  times = std::vector<double>();

  if (simparams.output_all_cycles) {
    int num_states =
        simparams.sim_num_time_steps / simparams.output_interval + 1;
    states.reserve(num_states);
    times.reserve(num_states);

  } else {
    int num_states =
        simparams.sim_pts_per_cycle / simparams.output_interval + 1;
    states.reserve(num_states);
    times.reserve(num_states);
  }
  double time = 0.0;

  // Run integrator
  DEBUG_MSG("Run time integration");
  int interval_counter = 0;
  int start_last_cycle =
      simparams.sim_num_time_steps - simparams.sim_pts_per_cycle;

  if (simparams.output_all_cycles || (0 >= start_last_cycle)) {
    times.push_back(time);
    states.push_back(std::move(state));
  }

  for (int i = 1; i < simparams.sim_num_time_steps; i++) {
    state = integrator.step(state, time);
    interval_counter += 1;
    time = simparams.sim_time_step_size * double(i);

    if ((interval_counter == simparams.output_interval) ||
        (!simparams.output_all_cycles && (i == start_last_cycle))) {
      if (simparams.output_all_cycles || (i >= start_last_cycle)) {
        times.push_back(time);
        states.push_back(std::move(state));
      }
      interval_counter = 0;
    }
  }

  DEBUG_MSG("Avg. number of nonlinear iterations per time step: " << integrator.avg_nonlin_iter());

  // Make times start from 0
  if (!simparams.output_all_cycles) {
    double start_time = times[0];
    for (auto& time : times) {
      time -= start_time;
    }
  }
}

std::vector<double> Solver::get_times() { return times; }

std::string Solver::get_full_result() {
  std::string output;

  if (simparams.output_variable_based) {
    output = to_variable_csv(times, states, model, simparams.output_mean_only,
                             simparams.output_derivative);

  } else {
    output = to_vessel_csv(times, states, model, simparams.output_mean_only,
                           simparams.output_derivative);
  }

  return output;
}

Eigen::VectorXd Solver::get_single_result(std::string dof_name) {
  int dof_index = model.dofhandler.get_variable_index(dof_name);
  int num_states = states.size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(num_states);

  for (size_t i = 0; i < num_states; i++) {
    result[i] = states[i].y[dof_index];
  }

  return result;
}

double Solver::get_single_result_avg(std::string dof_name) {
  int dof_index = model.dofhandler.get_variable_index(dof_name);
  int num_states = states.size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(num_states);

  for (size_t i = 0; i < num_states; i++) {
    result[i] = states[i].y[dof_index];
  }

  return result.mean();
}

void Solver::update_block_params(std::string block_name,
                                 std::vector<double> new_params) {
  auto block = model.get_block(block_name);

  if (new_params.size() != block->global_param_ids.size()) {
    throw std::runtime_error(
        "Parameter update failed! Number of provided parameters does not match "
        "with block parameters.");
  }

  for (size_t i = 0; i < new_params.size(); i++) {
    model.get_parameter(block->global_param_ids[i])->update(new_params[i]);
  }
}

void Solver::sanity_checks() {
  // Check that steady initial is not used with ClosedLoopHeartAndPulmonary
  if ((simparams.sim_steady_initial == true) &&
      (model.get_block("CLH") != nullptr)) {
    std::runtime_error(
        "ERROR: Steady initial condition is not compatible with "
        "ClosedLoopHeartAndPulmonary block.");
  }
}

void Solver::write_result_to_csv(std::string filename) {
  DEBUG_MSG("Write output");
  std::ofstream ofs(filename);
  ofs << get_full_result();
  ofs.close();
}
