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
 * @brief svZeroDSolver callinle interace.
 */
#include "interface.h"

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
//            Callible interface functions              //
//////////////////////////////////////////////////////////

extern "C" void initialize(const char* input_file, const double time_step, int& problem_id, int& system_size);

extern "C" void increment_time(const int problem_id, const double time, std::vector<double>& solution);

//------------
// initialize
//------------
// Initialize the 0D solver.
//
// Parameters:
//   input_file_arg: The name of the JSON 0D solver configuration file.
//   solver_time_step: The time step used by the 3D solver.
//
//   problem_id: The returned ID used to identify the 0D problem (block).
//
void initialize(const char* input_file_arg, const double solver_time_step, int& problem_id, int& system_size)
{
  std::cout << "========== svzero initialize ==========" << std::endl;
  std::string input_file(input_file_arg);
  std::cout << "[initialize] input_file: " << input_file << std::endl;
  std::string output_file = "svzerod.json";

  auto interface = new SolverInterface(input_file);
  interface->solver_time_step_ = solver_time_step;
  problem_id = interface->problem_id_;
  std::cout << "[initialize] problem_id: " << problem_id << std::endl;

  // Create configuration reader.
  IO::ConfigReader<T> config(input_file);

  // Create a model.
  auto model = config.get_model();
  interface->model_ = model; 
  system_size = model->dofhandler.size();
  std::cout << "[initialize] System size: " << system_size << std::endl;

  // Get simulation parameters from the input JSON file.
  //
  T time_step_size = config.get_time_step_size();
  int num_time_steps = config.get_num_time_steps();
  T absolute_tolerance = config.get_scalar_simulation_parameter("absolute_tolerance", 1e-8);
  int max_nliter = config.get_int_simulation_parameter("maximum_nonlinear_iterations", 30);
  int output_interval = config.get_int_simulation_parameter("output_interval", 1);
  bool steady_initial = config.get_bool_simulation_parameter("steady_initial", true);
  bool output_mean_only = config.get_bool_simulation_parameter("output_mean_only", false);

  std::cout << "[initialize] svzero parameters " << std::endl;
  std::cout << "[initialize]   num_time_steps: " << num_time_steps << std::endl;
  std::cout << "[initialize]   time_step_size: " << time_step_size << std::endl;
  std::cout << "[initialize]   absolute_tolerance: " << absolute_tolerance << std::endl;
  std::cout << "[initialize]   max_nliter: " << max_nliter << std::endl;

  interface->num_time_steps_ = num_time_steps;
  interface->time_step_size_ = time_step_size;
  interface->max_nliter_ = max_nliter;
  interface->absolute_tolerance_ = absolute_tolerance;

  // Create initial state. 
  ALGEBRA::State<T> state = ALGEBRA::State<T>::Zero(model->dofhandler.size());

  // Create steady initial state.
  //
  if (steady_initial) {
    std::cout << "[initialize] " << std::endl;
    std::cout << "[initialize] ----- Calculating steady initial condition ----- " << std::endl;
    T time_step_size_steady = config.cardiac_cycle_period / 10.0;
    std::cout << "[initialize] Create steady model ... " << std::endl;

    auto model_steady = config.get_model();
    model_steady->to_steady();

    ALGEBRA::Integrator<T, S> integrator_steady(*model_steady, time_step_size_steady, 0.1, absolute_tolerance, max_nliter);

    for (size_t i = 0; i < 31; i++) {
      state = integrator_steady.step(state, time_step_size_steady * T(i), *model_steady);
    }
  }

  interface->state_ = state;
  interface->time_step_ = 0;
  interface->save_interval_counter_ = 0;
  interface->output_interval_ = output_interval;

  std::cout << "[initialize] Done " << std::endl;
}

//----------------
// increment_time
//----------------
//
void increment_time(const int problem_id, const double time, std::vector<double>& solution)
{
  //std::cout << "[increment_time] " << std::endl;
  //std::cout << "[increment_time] ========== svzero increment_time ==========" << std::endl;
  auto interface = SolverInterface::interface_list_[problem_id];
  auto model = interface->model_;

  // [TODO] Use solver time step or svzero time step?
  auto time_step_size = interface->time_step_size_;
  auto absolute_tolerance = interface->absolute_tolerance_;
  auto max_nliter = interface->max_nliter_;
  //std::cout << "[increment_time] time: " << time << std::endl;

  ALGEBRA::Integrator<T, S> integrator(*model, time, 0.1, absolute_tolerance, max_nliter);
  auto state = interface->state_;
  interface->state_ = integrator.step(state, time, *model);
  interface->time_step_ += 1;

  int state_size = state.y.size();
  //std::cout << "[increment_time] state_size: " << state_size << std::endl;
  //std::cout << "[increment_time] solution size: " << solution.size() << std::endl;

  for (int i = 0; i < state.y.size(); i++) {
    solution[i] = state.y[i];
  }
}
