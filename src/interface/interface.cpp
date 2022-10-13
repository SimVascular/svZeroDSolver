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
 * @file interface.cpp
 * @brief svZeroDSolver callable interace.
 */
#include "interface.h"

#include <cmath>

typedef double T;

template <typename TT>
using S = ALGEBRA::SparseSystem<TT>;

// Static member data.
//
int SolverInterface::problem_id_count_ = 0;
std::map<int,SolverInterface*> SolverInterface::interface_list_;

//-----------------
// SolverInterface
//-----------------
//
SolverInterface::SolverInterface(const std::string& input_file_name) : input_file_name_(input_file_name)
{
  problem_id_ = problem_id_count_++;
  SolverInterface::interface_list_[problem_id_] = this;
}

SolverInterface::~SolverInterface()
{
}

//////////////////////////////////////////////////////////
//            Callable interface functions              //
//////////////////////////////////////////////////////////

extern "C" void initialize(std::string input_file, int& problem_id, int& pts_per_cycle, int& num_cycles, int& num_output_steps, std::vector<std::string>& block_names, std::vector<std::string>& variable_names);

extern "C" void increment_time(const int problem_id, const double external_time, std::vector<double>& solution);

extern "C" void run_simulation(const int problem_id, const double external_time, std::vector<double>& output_times, std::vector<double>& output_solutions, int& error_code);

extern "C" void update_block_params(const int problem_id, std::string block_name, std::vector<double>& params);

extern "C" void read_block_params(const int problem_id, std::string block_name, std::vector<double>& params);

extern "C" void get_block_node_IDs(const int problem_id, std::string block_name, std::vector<int>& IDs);

/**
 * @brief Initialize the 0D solver interface.
 *
 * @param input_file_arg The name of the JSON 0D solver configuration file.
 * @param external_time_step The time step used by the external program (3D solver).
 * @param problem_id The returned ID used to identify the 0D problem.
 * @param system_size Number of degrees-of-freedom.
 */
void initialize(std::string input_file_arg, int& problem_id, int& pts_per_cycle, int& num_cycles, int& num_output_steps, std::vector<std::string>& block_names, std::vector<std::string>& variable_names)
{
  DEBUG_MSG("========== svZeroD initialize ==========");
  std::cout << "[svZeroD::initialize] 1" << std::endl;
  std::string input_file(input_file_arg);
  DEBUG_MSG("[initialize] input_file: " << input_file);
  std::string output_file = "svzerod.csv";
  std::cout << "[svZeroD::initialize] 2" << std::endl;

  auto interface = new SolverInterface(input_file);
  problem_id = interface->problem_id_;
  DEBUG_MSG("[initialize] problem_id: " << problem_id);
  std::cout << "[svZeroD::initialize] 3" << std::endl;

  // Create configuration reader.
  IO::ConfigReader<T> reader;
  reader.load(input_file);

  // Create a model.
  auto model = reader.model;
  interface->model_ = model; 
  std::cout << "[svZeroD::initialize] 4" << std::endl;
//std::cout << "[initialize] block_index_map: " << model->block_index_map["branch0_seg0"] << std::endl;
//std::cout << "[initialize] block_index_map: " << model->block_index_map["INFLOW"] << std::endl;
//std::cout << "[initialize] block_index_map: " << model->block_index_map["OUT"] << std::endl;
  for(auto const &elem : model->block_index_map)
  {
    //std::cout << elem.first << " " << elem.second << "\n";
    block_names.push_back(elem.first);
  }
  std::cout << "[svZeroD::initialize] 5" << std::endl;
  variable_names = model->dofhandler.variables;
  std::cout << "[svZeroD::initialize] 6" << std::endl;

  // Get simulation parameters
  interface->time_step_size_ = reader.sim_time_step_size;
  interface->max_nliter_ = reader.sim_nliter;
  interface->absolute_tolerance_ = reader.sim_abs_tol;
  //interface->external_time_step_ = external_time_step;
  interface->time_step_ = 0;
  interface->system_size_ = model->dofhandler.size();
  //system_size = interface->system_size_;
  interface->output_interval_ = reader.output_interval;
  interface->num_time_steps_ = reader.sim_num_time_steps;
  interface->pts_per_cycle_ = reader.sim_pts_per_cycle;
  pts_per_cycle = reader.sim_pts_per_cycle;
  num_cycles = reader.sim_num_cycles;

  std::cout << "[svZeroD::initialize] 7" << std::endl;
  // For how many time steps are outputs being returned?
  if (reader.output_mean_only) {
    num_output_steps = 1;
  } else if (reader.output_last_cycle_only) {
    num_output_steps = interface->pts_per_cycle_;
  } else {
    //num_output_steps = std::ceil(interface->num_time_steps_/interface->output_interval_);
    num_output_steps = interface->num_time_steps_;
  }
  interface->num_output_steps_ = num_output_steps;
  DEBUG_MSG("[initialize] System size: " << interface->system_size_);

  // Create initial state.
  ALGEBRA::State<T> state = reader.initial_state;

  // Create steady initial state.
  if (reader.sim_steady_initial) {
    DEBUG_MSG("[initialize] ----- Calculating steady initial condition ----- ");
    T time_step_size_steady = reader.sim_cardiac_cycle_period / 10.0;
    DEBUG_MSG("[initialize] Create steady model ... ");

    auto model_steady = reader.model;
    model_steady->to_steady();

    ALGEBRA::Integrator<T> integrator_steady(*model_steady, time_step_size_steady, 0.1, interface->absolute_tolerance_, interface->max_nliter_);

    for (size_t i = 0; i < 31; i++) {
      state = integrator_steady.step(state, time_step_size_steady * T(i), *model_steady);
    }
  }
  interface->state_ = state;

  DEBUG_MSG("[initialize] Done");
  std::cout << "[svZeroD::initialize] END" << std::endl;
}

/**
 * @brief Update the parameters of a particular block.
 *
 * @param problem_id The returned ID used to identify the 0D problem.
 * @param block name The name of the block to update.
 * @param params New parameters for the block (structure depends on block type).
 */
void update_block_params(const int problem_id, std::string block_name, std::vector<double>& params)
{
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;
  std::cout << "[update_block_params] block_index_map: " << model->block_index_map["branch0_seg0"] << std::endl;
  std::cout << "[update_block_params] block_index_map: " << model->block_index_map["INFLOW"] << std::endl;
  std::cout << "[update_block_params] block_index_map: " << model->block_index_map["OUT"] << std::endl;
  for(auto const &elem : model->block_index_map)
  {
    std::cout << elem.first << " " << elem.second << "\n";
  }
  int block_index = model->block_index_map[block_name];
  std::cout << "[update_block_params] Input block name: " << block_name << std::endl;
  std::cout << "[update_block_params] block_index: " << block_index << std::endl;
  auto block = model->blocks[block_index];
  std::cout << "[update_block_params] block name: " << block->name << std::endl;
  block->update_block_params(params);
}

/**
 * @brief Read the parameters of a particular block.
 *
 * @param problem_id The returned ID used to identify the 0D problem.
 * @param block name The name of the block to read.
 * @param params Parameters of the block (structure depends on block type).
 */
//void read_block_params(const int problem_id, const char* block_name, std::vector<double>& params)
void read_block_params(const int problem_id, std::string block_name, std::vector<double>& params)
{
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;
  int block_index = model->block_index_map[block_name];
  auto block = model->blocks[block_index];
  block->get_block_params(params);
}

/**
 * @brief Return the IDs of the input and output node in the solution vector.
 *
 * @param problem_id The returned ID used to identify the 0D problem.
 * @param block name The name of the block whose node IDs are returned.
 * @param IDs Vector containing IDs of input and output nodes in the following order: {num inlet nodes, inlet flow[0], inlet pressure[0],..., num outlet nodes, outlet flow[0], outlet pressure[0],...}.
 * @param specifier String containing any extra information. For external coupling blocks, this stores location of the block (inlet/outlet).
 */
void get_block_node_IDs(const int problem_id, std::string block_name, std::vector<int>& IDs)
{
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;
  int block_index = model->block_index_map[block_name];
  auto block = model->blocks[block_index];
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
 * @brief Increment the 0D solution by one time step.
 *
 * @param problem_id The ID used to identify the 0D problem.
 * @param external_time The current time in the external program.
 * @param solution The solution vector containing all degrees-of-freedom.
 */
void increment_time(const int problem_id, const double external_time, std::vector<double>& solution)
{
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;

  auto time_step_size = interface->time_step_size_;
  auto absolute_tolerance = interface->absolute_tolerance_;
  auto max_nliter = interface->max_nliter_;

  ALGEBRA::Integrator<T> integrator(*model, time_step_size, 0.1, absolute_tolerance, max_nliter);
  auto state = interface->state_;
  interface->state_ = integrator.step(state, external_time, *model);
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
 * @param output_times Vector containing time-stamps for output_solution vectors.
 * @param output_solutions The solution vector containing all degrees-of-freedom stored sequentially (1D vector).
 */
void run_simulation(const int problem_id, const double external_time, std::vector<double>& output_times, std::vector<double>& output_solutions, int& error_code)
{
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;
//std::cout << "[run_simulation] block_index_map: " << model->block_index_map["branch0_seg0"] << std::endl;
//std::cout << "[run_simulation] block_index_map: " << model->block_index_map["INFLOW"] << std::endl;
//std::cout << "[run_simulation] block_index_map: " << model->block_index_map["OUT"] << std::endl;

  auto time_step_size = interface->time_step_size_;
  auto absolute_tolerance = interface->absolute_tolerance_;
  auto max_nliter = interface->max_nliter_;
  auto num_time_steps = interface->num_time_steps_;
  auto system_size = interface->system_size_;
  auto num_output_steps = interface->num_output_steps_;

  //std::cout << "[run_simulation] 1 " << std::endl;

  ALGEBRA::Integrator<T> integrator(*model, time_step_size, 0.1, absolute_tolerance, max_nliter);
  auto state = interface->state_;
  T time = external_time;
  std::vector<T> times;
  times.reserve(num_time_steps);
  times.push_back(time);
  std::vector<ALGEBRA::State<T>> states;
  states.reserve(num_time_steps);
  states.push_back(state);

  //std::cout << "[run_simulation] 2 " << std::endl;
  // Run integrator
  interface->time_step_ = 0;
  error_code = 0;
  bool isNaN = false;
  for (int i = 1; i < num_time_steps; i++) {
    //std::cout << "[run_simulation] time: " << time << std::endl;
    interface->time_step_ += 1;
    state = integrator.step(state, time, *model);
    // Check for NaNs in the state vector
    if ((i % 100) == 0) {
      //std::cout << "[run_simulation] time, isNaN: " << time <<", "<<isNaN<< std::endl;
      for(int j=0; j < system_size; j++) {
        isNaN = (state.y[j] != state.y[j]);
        if(isNaN) {
          std::cout << "Found NaN in state vector at timestep " << i << " and index " << j << std::endl;
          error_code = 1;
          return;
        }
      }
    }
    time += time_step_size;
    times.push_back(time);
    states.push_back(std::move(state));
  }
  //std::cout << "[run_simulation] 3 " << std::endl;

  // Extract last cardiac cycle
  if (interface->output_last_cycle_only_) {
    states.erase(states.begin(), states.end() - interface->pts_per_cycle_);
    times.erase(times.begin(), times.end() - interface->pts_per_cycle_);
    T start_time = times[0];
    for (auto& time : times) {
      time -= start_time;
    }
  }

  // Write states to solution output vector
  if (output_solutions.size() != num_output_steps*system_size) {
    throw std::runtime_error("Solution vector size is wrong.");
  }
  int output_idx = 0;
  int soln_idx = 0;
  for (int t = 0; t < num_output_steps; t++) {
    auto state = states[t];
    output_times[t] = times[t];
    for (int i = 0; i < system_size; i++) {
      soln_idx = output_idx*system_size + i;
      output_solutions[soln_idx] = state.y[i];
    }
    output_idx++;
  }
  
//output_idx = 0;
//soln_idx = 0;
//for (int t = 0; t < num_output_steps; t++) {
//  std::cout << "[run_simulation] output_times: " << output_times[t] << " output_solutions: ";
//  for (int i = 0; i < system_size; i++) {
//    soln_idx = system_size*output_idx + i;
//    std::cout << output_solutions[soln_idx] << " ";
//  }
//  std::cout << std::endl;
//  output_idx++;
//}
}
