#include "Solver.h"

#include "csv_writer.h"

Solver::Solver(const nlohmann::json& config) {
  validate_input(config);
  DEBUG_MSG("Read simulation parameters");
  simparams = load_simulation_params(config);
  DEBUG_MSG("Load model");
  //model = Model();
  this->model = std::shared_ptr<Model>(new Model());
  //load_simulation_model(config, model);
  load_simulation_model(config, *this->model.get());
  DEBUG_MSG("Load initial condition");
  //initial_state = load_initial_condition(config, model);
  initial_state = load_initial_condition(config, *this->model.get());

  //DEBUG_MSG("Cardiac cycle period " << model.cardiac_cycle_period);
  DEBUG_MSG("Cardiac cycle period " << this->model->cardiac_cycle_period);

  // Calculate time step size
  if (!simparams.sim_coupled) {
    //simparams.sim_time_step_size = model.cardiac_cycle_period /
    simparams.sim_time_step_size = this->model->cardiac_cycle_period /
                                   (double(simparams.sim_pts_per_cycle) - 1.0);
  } else {
    simparams.sim_time_step_size = simparams.sim_external_step_size /
                                   (double(simparams.sim_num_time_steps) - 1.0);
  }
  //std::cout<<"[Solver] cardiac_cycle_period = "<<model.cardiac_cycle_period<<std::endl;
  std::cout<<"[Solver] cardiac_cycle_period = "<<this->model->cardiac_cycle_period<<std::endl;
  std::cout<<"[Solver] test0"<<std::endl;
  sanity_checks();
  std::cout<<"[Solver] test1"<<std::endl;
  //this->model = model;
  std::cout<<"[Solver] test2"<<std::endl;
}

// Solver::~Solver() {}

void Solver::run() {
  auto state = initial_state;

  std::cout<<"test1"<<std::endl;
  std::cout<<simparams.sim_steady_initial<<std::endl;
  //std::cout<<"[run] cardiac_cycle_period = "<<model.cardiac_cycle_period<<std::endl;
  std::cout<<"[run] cardiac_cycle_period = "<<this->model->cardiac_cycle_period<<std::endl;
  // Create steady initial
  if (simparams.sim_steady_initial) {
    DEBUG_MSG("Calculate steady initial condition");
    //double time_step_size_steady = model.cardiac_cycle_period / 10.0;
    double time_step_size_steady = this->model->cardiac_cycle_period / 10.0;
    //auto model_steady = this->model;
    //model_steady->to_steady();
    this->model->to_steady();

//  Integrator integrator_steady(&model, time_step_size_steady,
//                               simparams.sim_rho_infty, simparams.sim_abs_tol,
//                               simparams.sim_nliter);
    Integrator integrator_steady(this->model.get(), time_step_size_steady,
                                 simparams.sim_rho_infty, simparams.sim_abs_tol,
                                 simparams.sim_nliter);

    for (int i = 0; i < 31; i++) {
      state = integrator_steady.step(state, time_step_size_steady * double(i));
    }

    //model.to_unsteady();
    this->model->to_unsteady();
  }
  std::cout<<"test2"<<std::endl;

  // Set-up integrator
  DEBUG_MSG("Setup time integration");
//Integrator integrator(&model, simparams.sim_time_step_size,
//                      simparams.sim_rho_infty, simparams.sim_abs_tol,
//                      simparams.sim_nliter);
  Integrator integrator(this->model.get(), simparams.sim_time_step_size,
                        simparams.sim_rho_infty, simparams.sim_abs_tol,
                        simparams.sim_nliter);
  std::cout<<"test3"<<std::endl;

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
  std::cout<<"test4"<<std::endl;
  std::cout<<simparams.output_all_cycles<<" "<<simparams.sim_num_time_steps<<" "<<simparams.output_interval<<" "<<simparams.sim_pts_per_cycle<<std::endl;

  // Run integrator
  DEBUG_MSG("Run time integration");
  int interval_counter = 0;
  int start_last_cycle =
      simparams.sim_num_time_steps - simparams.sim_pts_per_cycle;

  if (simparams.output_all_cycles || (0 >= start_last_cycle)) {
    times.push_back(time);
    states.push_back(std::move(state));
  }
  std::cout<<"test5"<<std::endl;
  std::cout<<time<<std::endl;
//for (int i = 0; i<23; i++) {
//  std::cout<<state.y[i]<<" ";
//}
  std::cout<<std::endl;

  for (int i = 1; i < simparams.sim_num_time_steps; i++) {
    if ( i < 2) {
        std::cout<<i<<" ";
//      for (int j = 0; j<23; j++) {
//        std::cout<<state.y[j]<<" ";
//      }
        std::cout<<std::endl;
    }
    //std::cout<<"test"<<" ";
    //std::cout<<time<<" ";
    state = integrator.step(state, time);
    //std::cout<<"test-after step"<<std::endl;
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
    //std::cout<<"test-loop"<<std::endl;
  }
  std::cout<<std::endl;

  DEBUG_MSG("Avg. number of nonlinear iterations per time step: "
            << integrator.avg_nonlin_iter());
  std::cout<<"test6"<<std::endl;

  // Make times start from 0
  if (!simparams.output_all_cycles) {
    double start_time = times[0];
    for (auto& time : times) {
      time -= start_time;
    }
  }
}

std::vector<double> Solver::get_times() const { return times; }

std::string Solver::get_full_result() const {
  std::string output;

  if (simparams.output_variable_based) {
//  output = to_variable_csv(times, states, model, simparams.output_mean_only,
//                           simparams.output_derivative);
    output = to_variable_csv(times, states, *this->model.get(), simparams.output_mean_only,
                             simparams.output_derivative);


  } else {
//  output = to_vessel_csv(times, states, model, simparams.output_mean_only,
//                         simparams.output_derivative);
    output = to_vessel_csv(times, states, *this->model.get(), simparams.output_mean_only,
                           simparams.output_derivative);
  }

  return output;
}

Eigen::VectorXd Solver::get_single_result(const std::string& dof_name) const {
  //int dof_index = model.dofhandler.get_variable_index(dof_name);
  int dof_index = this->model->dofhandler.get_variable_index(dof_name);
  int num_states = states.size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(num_states);

  for (size_t i = 0; i < num_states; i++) {
    result[i] = states[i].y[dof_index];
  }

  return result;
}

double Solver::get_single_result_avg(const std::string& dof_name) const {
  //int dof_index = model.dofhandler.get_variable_index(dof_name);
  int dof_index = this->model->dofhandler.get_variable_index(dof_name);
  int num_states = states.size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(num_states);

  for (size_t i = 0; i < num_states; i++) {
    result[i] = states[i].y[dof_index];
  }

  return result.mean();
}

void Solver::update_block_params(const std::string& block_name,
                                 const std::vector<double>& new_params) {
  //auto block = model.get_block(block_name);
  auto block = this->model->get_block(block_name);

  if (new_params.size() != block->global_param_ids.size()) {
    throw std::runtime_error(
        "Parameter update failed! Number of provided parameters does not match "
        "with block parameters.");
  }

  for (size_t i = 0; i < new_params.size(); i++) {
    //model.get_parameter(block->global_param_ids[i])->update(new_params[i]);
    this->model->get_parameter(block->global_param_ids[i])->update(new_params[i]);
    // parameter_values vector needs to be seperately updated for constant parameters. This does not need to be done for time-dependent parameters because it is handled in Model::update_time
    //model.update_parameter_value(block->global_param_ids[i], new_params[i]);
    this->model->update_parameter_value(block->global_param_ids[i], new_params[i]);
  }
}

void Solver::sanity_checks() {
  std::cout<<"[sanity_check] test0"<<std::endl;
  // Check that steady initial is not used with ClosedLoopHeartAndPulmonary
  //if ((simparams.sim_steady_initial == true) && (model.has_block("CLH"))) {
  if ((simparams.sim_steady_initial == true) && (this->model->has_block("CLH"))) {
    std::runtime_error(
        "ERROR: Steady initial condition is not compatible with "
        "ClosedLoopHeartAndPulmonary block.");
  }
  std::cout<<"[sanity_check] test1"<<std::endl;
}

void Solver::write_result_to_csv(const std::string& filename) const {
  DEBUG_MSG("Write output");
  std::ofstream ofs(filename);
  ofs << get_full_result();
  ofs.close();
}
