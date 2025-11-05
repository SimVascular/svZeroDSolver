// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "interface.h"

#include <cmath>

#if defined(_WIN32) || defined(__CYGWIN__)
#define SVZEROD_INTERFACE_API __declspec(dllexport)
#else
#define SVZEROD_INTERFACE_API
#endif

#include "SimulationParameters.h"

// Static member data.
int SolverInterface::problem_id_count_ = 0;
std::map<int, SolverInterface*> SolverInterface::interface_list_;

//-----------------
// SolverInterface
//-----------------
SolverInterface::SolverInterface(const std::string& input_file_name)
    : input_file_name_(input_file_name) {
  problem_id_ = problem_id_count_++;
  SolverInterface::interface_list_[problem_id_] = this;
}

SolverInterface::~SolverInterface() {}

//////////////////////////////////////////////////////////
//            Callable interface functions              //
//////////////////////////////////////////////////////////

/**
 * @brief Initialize the 0D solver interface.
 *
 * @param input_file_arg The name of the JSON 0D solver configuration file.
 * @param problem_id The returned ID used to identify the 0D problem.
 * @param pts_per_cycle Number of time steps per cycle in the 0D model.
 * @param num_cycles Number of cardiac cycles in the 0D model.
 * @param num_output_steps Number of steps at which outputs are recorded.
 * @param block_names Vector of all the 0D block names.
 * @param variable_names Vector of all the 0D variable names.
 */
extern "C" SVZEROD_INTERFACE_API void initialize(
    const char* input_file_arg, int& problem_id, int& pts_per_cycle,
    int& num_cycles, int& num_output_steps,
    std::vector<std::string>& block_names,
    std::vector<std::string>& variable_names) {
  // Print immediately to see if we enter the function
  std::cerr << "[DLL:initialize] ENTRY POINT REACHED\n";
  std::cerr.flush();

  DEBUG_MSG("========== svZeroD initialize ==========");

  // WORKAROUND: Copy strings element by element using C strings to avoid vector reallocation
  // We'll populate the vectors without calling resize/clear/etc
  std::cerr << "[DLL:initialize] Preparing to populate vectors safely\n";
  std::cerr.flush();

  // Convert C string to std::string inside DLL to avoid ABI issues
  std::string input_file(input_file_arg ? input_file_arg : "");
  DEBUG_MSG("[initialize] input_file: " << input_file);

  // Always print this to help debug Release builds
  std::cerr << "[DLL:initialize] Received file path: " << input_file << std::endl;
  std::cerr << "[DLL:initialize] Path length: " << input_file.length() << std::endl;
  std::cerr.flush();

  std::string output_file = "svzerod.csv";

  auto interface = new SolverInterface(input_file);
  problem_id = interface->problem_id_;
  DEBUG_MSG("[initialize] problem_id: " << problem_id);

  // Create configuration reader.
  std::cerr << "[DLL:initialize] Opening file..." << std::endl;
  std::cerr.flush();
  std::ifstream ifs(input_file);
  if (!ifs.is_open() || !ifs.good()) {
    std::cerr << "[DLL:initialize] FAILED to open file!" << std::endl;
    std::cerr.flush();
    throw std::runtime_error("Failed to open input file: " + input_file);
  }
  std::cerr << "[DLL:initialize] File opened successfully" << std::endl;
  std::cerr.flush();
  const auto& config = nlohmann::json::parse(ifs);
  auto simparams = load_simulation_params(config);

  auto model = std::shared_ptr<Model>(new Model());

  load_simulation_model(config, *model.get());
  auto state = load_initial_condition(config, *model.get());

  // Check that steady initial is not set when ClosedLoopHeartAndPulmonary is
  // used
  if ((simparams.sim_steady_initial == true) && (model->has_block("CLH"))) {
    std::runtime_error(
        "ERROR: Steady initial condition is not compatible with "
        "ClosedLoopHeartAndPulmonary block.");
  }

  // Set default cardiac cycle period if not set by model
  if (model->cardiac_cycle_period < 0.0) {
    model->cardiac_cycle_period =
        1.0;  // If it has not been read from config or Parameter
              // yet, set as default value of 1.0
  }

  // Calculate time step size
  if (!simparams.sim_coupled) {
    simparams.sim_time_step_size = model->cardiac_cycle_period /
                                   (double(simparams.sim_pts_per_cycle) - 1.0);
  } else {
    simparams.sim_time_step_size = simparams.sim_external_step_size /
                                   (double(simparams.sim_num_time_steps) - 1.0);
  }

  // Create a model.
  interface->model_ = model;

  std::cerr << "[DLL:initialize] About to populate internal vectors\n";
  std::cerr.flush();

  // Populate DLL-internal vectors (safe because they're allocated in DLL)
  interface->block_names_.clear();
  interface->variable_names_.clear();

  for (size_t i = 0; i < model->get_num_blocks(); i++) {
    interface->block_names_.push_back(model->get_block(i)->get_name());
  }
  interface->variable_names_ = model->dofhandler.variables;

  std::cerr << "[DLL:initialize] Internal vectors populated, copying to caller\n";
  std::cerr << "[DLL:initialize] block_names size: " << interface->block_names_.size() << "\n";
  std::cerr << "[DLL:initialize] variable_names size: " << interface->variable_names_.size() << "\n";
  std::cerr.flush();

  // Copy to caller's vectors element-by-element using indexed assignment (safest way)
  // This avoids any resizing operations on the caller's vectors
  for (size_t i = 0; i < interface->block_names_.size(); i++) {
    if (i < block_names.size()) {
      // Assign to existing element
      block_names[i] = interface->block_names_[i].c_str();
    } else {
      // Need to grow - try push_back
      try {
        block_names.push_back(interface->block_names_[i].c_str());
      } catch (...) {
        std::cerr << "[DLL:initialize] Failed to push_back block_name " << i << "\n";
        std::cerr.flush();
        // Can't modify caller's vector - just skip
        break;
      }
    }
  }

  for (size_t i = 0; i < interface->variable_names_.size(); i++) {
    if (i < variable_names.size()) {
      variable_names[i] = interface->variable_names_[i].c_str();
    } else {
      try {
        variable_names.push_back(interface->variable_names_[i].c_str());
      } catch (...) {
        std::cerr << "[DLL:initialize] Failed to push_back variable_name " << i << "\n";
        std::cerr.flush();
        break;
      }
    }
  }

  std::cerr << "[DLL:initialize] Finished copying to caller's vectors\n";
  std::cerr.flush();

  // Get simulation parameters
  interface->time_step_size_ = simparams.sim_time_step_size;
  interface->rho_infty_ = simparams.sim_rho_infty;
  interface->max_nliter_ = simparams.sim_nliter;
  interface->absolute_tolerance_ = simparams.sim_abs_tol;
  interface->time_step_ = 0;
  interface->system_size_ = model->dofhandler.size();
  interface->num_time_steps_ = simparams.sim_num_time_steps;
  interface->pts_per_cycle_ = simparams.sim_pts_per_cycle;
  pts_per_cycle = simparams.sim_pts_per_cycle;
  num_cycles = simparams.sim_num_cycles;
  interface->external_step_size_ = simparams.sim_external_step_size;

  // For how many time steps are outputs being returned?
  // NOTE: Only tested num_output_steps = interface->num_time_steps_
  if (simparams.output_mean_only) {
    num_output_steps = 1;
    throw std::runtime_error(
        "ERROR: Option output_mean_only has not been implemented when "
        "using the svZeroDSolver interface library.");
  } else if (!simparams.output_all_cycles) {
    num_output_steps = interface->pts_per_cycle_;
    throw std::runtime_error(
        "ERROR: Option output_last_cycle_only has been implemented but not "
        "tested when using the svZeroDSolver interface library. Please test "
        "this "
        "functionality before removing this message.");
  } else {
    num_output_steps = interface->num_time_steps_;
  }
  interface->num_output_steps_ = num_output_steps;
  DEBUG_MSG("[initialize] System size: " << interface->system_size_);

  // Create steady initial state.
  if (simparams.sim_steady_initial) {
    DEBUG_MSG("[initialize] ----- Calculating steady initial condition ----- ");
    double time_step_size_steady = model->cardiac_cycle_period / 10.0;
    DEBUG_MSG("[initialize] Create steady model ... ");

    auto model_steady = model;
    model_steady->to_steady();
    Integrator integrator_steady(
        model_steady.get(), time_step_size_steady, interface->rho_infty_,
        interface->absolute_tolerance_, interface->max_nliter_);

    for (size_t i = 0; i < 31; i++) {
      state = integrator_steady.step(state, time_step_size_steady * double(i));
    }
    model_steady->to_unsteady();
  }

  interface->state_ = state;

  // Initialize states and times vectors because size is now known
  interface->times_.resize(num_output_steps);
  interface->states_.resize(num_output_steps);

  // Initialize integrator
  interface->integrator_ =
      Integrator(model.get(), interface->time_step_size_, interface->rho_infty_,
                 interface->absolute_tolerance_, interface->max_nliter_);

  DEBUG_MSG("[initialize] Done");
}

/**
 * @brief Set the timestep of the external program. For cases when 0D time step
 * depends on external time step.
 *
 * @param problem_id The returned ID used to identify the 0D problem.
 * @param external_step_size The time step size of the external program.
 */
extern "C" SVZEROD_INTERFACE_API void set_external_step_size(
    int problem_id, double external_step_size) {
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;

  // Update external step size in model and interface
  interface->external_step_size_ = external_step_size;

  // Update time step size in interface
  double zerod_step_size =
      external_step_size / (double(interface->num_time_steps_) - 1.0);
  interface->time_step_size_ = zerod_step_size;
}

/**
 * @brief Update the parameters of a particular block.
 *
 * @param problem_id The returned ID used to identify the 0D problem.
 * @param block name The name of the block to update.
 * @param params New parameters for the block (structure depends on block type).
 */
extern "C" SVZEROD_INTERFACE_API void update_block_params(
    int problem_id, const char* block_name_arg,
    std::vector<double>& params) {
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;

  // Convert C string to std::string inside DLL to avoid ABI issues
  std::string block_name(block_name_arg ? block_name_arg : "");

  // Find the required block
  auto block = model->get_block(block_name);
  auto block_type = model->get_block_type(block_name);
  // Update is handled differently for blocks that have time-varying parameters
  // (PRESSUREBC and FLOWBC)
  // TODO: Does this need to be done for OPENLOOPCORONARYBC and RESISTANCEBC
  // too?
  if ((block_type == BlockType::pressure_bc) ||
      (block_type == BlockType::flow_bc)) {
    std::vector<double> times_new;
    std::vector<double> values_new;
    int num_time_pts = (int)params[0];
    for (int i = 0; i < num_time_pts; i++) {
      times_new.push_back(params[1 + i]);
      values_new.push_back(params[1 + num_time_pts + i]);
    }
    model->get_parameter(block->global_param_ids[0])
        ->update(times_new, values_new);
  } else {
    if (block->global_param_ids.size() != params.size()) {
      throw std::runtime_error(
          "New parameter vector (given size = " +
          std::to_string(params.size()) +
          ") does not match number of parameters of block " + block_name +
          " (required size = " +
          std::to_string(block->global_param_ids.size()) + ")");
    }
    for (size_t i = 0; i < params.size(); i++) {
      model->get_parameter(block->global_param_ids[i])->update(params[i]);
      // parameter_values vector needs to be seperately updated for constant
      // parameters. This does not need to be done for time-dependent parameters
      // because it is handled in Model::update_time
      model->update_parameter_value(block->global_param_ids[i], params[i]);
    }
  }
}

/**
 * @brief Read the parameters of a particular block.
 *
 * @param problem_id The returned ID used to identify the 0D problem.
 * @param block name The name of the block to read.
 * @param params Parameters of the block (structure depends on block type).
 */
extern "C" SVZEROD_INTERFACE_API void read_block_params(
    int problem_id, const char* block_name_arg,
    std::vector<double>& params) {
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;

  // Convert C string to std::string inside DLL to avoid ABI issues
  std::string block_name(block_name_arg ? block_name_arg : "");

  auto block = model->get_block(block_name);
  if (params.size() != block->global_param_ids.size()) {
    throw std::runtime_error(
        "Parameter vector (given size = " + std::to_string(params.size()) +
        ") does not match number of parameters of block " + block_name +
        " (required size = " + std::to_string(block->global_param_ids.size()) +
        ")");
  }
  for (size_t i = 0; i < params.size(); i++) {
    params[i] = model->get_parameter_value(block->global_param_ids[i]);
  }
}

/**
 * @brief Return the IDs of the input and output nodes in the solution vector
 * for a given block.
 *
 * @param problem_id The returned ID used to identify the 0D problem.
 * @param block name The name of the block whose node IDs are returned.
 * @param IDs Vector containing IDs of input and output nodes in the following
 * order: {num inlet nodes, inlet flow[0], inlet pressure[0],..., num outlet
 * nodes, outlet flow[0], outlet pressure[0],...}.
 */
extern "C" SVZEROD_INTERFACE_API void get_block_node_IDs(
    int problem_id, const char* block_name_arg, std::vector<int>& IDs) {
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;

  // Convert C string to std::string inside DLL to avoid ABI issues
  std::string block_name(block_name_arg ? block_name_arg : "");

  // Find the required block
  auto block = model->get_block(block_name);

  // IDs are stored in the following format
  // {num inlet nodes, inlet flow[0], inlet pressure[0],..., num outlet nodes,
  // outlet flow[0], outlet pressure[0],...}
  IDs.clear();
  IDs.push_back(block->inlet_nodes.size());
  for (int i = 0; i < block->inlet_nodes.size(); i++) {
    IDs.push_back(block->inlet_nodes[i]->flow_dof);
    IDs.push_back(block->inlet_nodes[i]->pres_dof);
  }
  IDs.push_back(block->outlet_nodes.size());
  for (int i = 0; i < block->outlet_nodes.size(); i++) {
    IDs.push_back(block->outlet_nodes[i]->flow_dof);
    IDs.push_back(block->outlet_nodes[i]->pres_dof);
  }
}

/**
 * @brief Return the y state vector.
 *
 * @param problem_id The ID used to identify the 0D problem.
 * @param y The state vector containing all state.y degrees-of-freedom.
 */
extern "C" SVZEROD_INTERFACE_API void return_y(int problem_id,
                                               std::vector<double>& y) {
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;
  auto system_size = interface->system_size_;
  if (y.size() != system_size) {
    throw std::runtime_error(
        "ERROR: State vector size is wrong in return_y().");
  }

  auto state = interface->state_;
  for (int i = 0; i < system_size; i++) {
    y[i] = state.y[i];
  }
}

/**
 * @brief Return the ydot state vector.
 *
 * @param problem_id The ID used to identify the 0D problem.
 * @param ydot The state vector containing all state.ydot degrees-of-freedom.
 */
extern "C" SVZEROD_INTERFACE_API void return_ydot(int problem_id,
                                                  std::vector<double>& ydot) {
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;
  auto system_size = interface->system_size_;
  if (ydot.size() != system_size) {
    throw std::runtime_error(
        "ERROR: State vector size is wrong in return_ydot().");
  }

  auto state = interface->state_;
  for (int i = 0; i < system_size; i++) {
    ydot[i] = state.ydot[i];
  }
}

/**
 * @brief Update the state vector.
 *
 * @param problem_id The ID used to identify the 0D problem.
 * @param new_state_y The new state vector containing all state.y
 * degrees-of-freedom.
 * @param new_state_ydot The new state vector containing all state.ydot
 * degrees-of-freedom.
 */
extern "C" SVZEROD_INTERFACE_API void update_state(
    int problem_id, const std::vector<double>& new_state_y,
    const std::vector<double>& new_state_ydot) {
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;
  auto system_size = interface->system_size_;
  if ((new_state_y.size() != system_size) ||
      (new_state_ydot.size() != system_size)) {
    throw std::runtime_error(
        "ERROR: State vector size is wrong in update_state().");
  }

  auto state = interface->state_;
  for (int i = 0; i < system_size; i++) {
    state.y[i] = new_state_y[i];
    state.ydot[i] = new_state_ydot[i];
  }
  interface->state_ = state;
}

/**
 * @brief Increment the 0D solution by one time step.
 *
 * @param problem_id The ID used to identify the 0D problem.
 * @param external_time The current time in the external program.
 * @param solution The solution vector containing all degrees-of-freedom.
 */
extern "C" SVZEROD_INTERFACE_API void increment_time(
    int problem_id, const double external_time, std::vector<double>& solution) {
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;

  auto time_step_size = interface->time_step_size_;
  auto absolute_tolerance = interface->absolute_tolerance_;
  auto max_nliter = interface->max_nliter_;
  Integrator integrator(model.get(), time_step_size, interface->rho_infty_,
                        absolute_tolerance, max_nliter);
  auto state = interface->state_;
  interface->state_ = integrator.step(state, external_time);
  interface->time_step_ += 1;

  for (int i = 0; i < state.y.size(); i++) {
    solution[i] = state.y[i];
  }
}

/**
 * @brief Run a full 0D simululation.
 *
 * @param problem_id The ID used to identify the 0D problem.
 * @param external_time The current time in the external program.
 * @param output_times Vector containing time-stamps for output_solution
 * vectors.
 * @param output_solutions The solution vector containing all degrees-of-freedom
 * stored sequentially (1D vector).
 * @param error_code This is 1 if a NaN is found in the solution vector, 0
 * otherwise.
 */
extern "C" SVZEROD_INTERFACE_API void run_simulation(
    int problem_id, const double external_time,
    std::vector<double>& output_times, std::vector<double>& output_solutions,
    int& error_code) {
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;

  auto time_step_size = interface->time_step_size_;
  auto absolute_tolerance = interface->absolute_tolerance_;
  auto max_nliter = interface->max_nliter_;
  auto num_time_steps = interface->num_time_steps_;
  auto system_size = interface->system_size_;
  auto num_output_steps = interface->num_output_steps_;

  auto integrator = interface->integrator_;
  integrator.update_params(time_step_size);

  auto state = interface->state_;
  double time = external_time;

  interface->times_[0] = time;
  interface->states_[0] = state;

  // Run integrator
  interface->time_step_ = 0;
  error_code = 0;
  bool isNaN = false;
  for (int i = 1; i < num_time_steps; i++) {
    interface->time_step_ += 1;
    state = integrator.step(state, time);
    // Check for NaNs in the state vector
    if ((i % 100) == 0) {
      for (int j = 0; j < system_size; j++) {
        isNaN = (state.y[j] != state.y[j]);
        if (isNaN) {
          std::cout << "Found NaN in state vector at timestep " << i
                    << " and index " << j << std::endl;
          error_code = 1;
          return;
        }
      }
    }
    time += time_step_size;
    interface->times_[i] = time;
    interface->states_[i] = state;
  }
  interface->state_ = state;

  // Write states to solution output vector
  if (output_solutions.size() != num_output_steps * system_size) {
    throw std::runtime_error("Solution vector size is wrong.");
  }
  int output_idx = 0;
  int soln_idx = 0;
  int start_idx = 0;
  double start_time = 0.0;
  if (interface->output_last_cycle_only_) {  // NOT TESTED
    start_idx = interface->num_time_steps_ - interface->pts_per_cycle_;
    start_time = interface->times_[start_idx];
    throw std::runtime_error(
        "ERROR: Option output_last_cycle_only has been implemented but not "
        "tested when using the svZeroDSolver interface library. Please test "
        "this "
        "functionality before removing this message.");
  }
  for (int t = start_idx; t < num_output_steps; t++) {
    state = interface->states_[t];
    output_times[t] = interface->times_[t] - start_time;
    for (int i = 0; i < system_size; i++) {
      soln_idx = output_idx * system_size + i;
      output_solutions[soln_idx] = state.y[i];
    }
    output_idx++;
  }

  // Release dynamic memory
  // integrator.clean();
}
