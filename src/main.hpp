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
// clang-format off
/**
 * @file main.hpp
 * @brief Main page documentation
 *
 * @mainpage svZeroDSolver
 *
 * \tableofcontents
 *
 * The svZeroDSolver is a fast simulation tool for modeling the hemodynamics of
 * vascular networks using zero-dimensional (0D) lumped parameter models.
 *
 * * <a href="https://github.com/SimVascular/svZeroDSolver">Source
 * repository</a>
 * * <a href="https://simvascular.github.io">About SimVascular</a>
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R$, *-] (3,0)
 * to [L, l=$L$, *-*] (5,0)
 * node[anchor=south]{$P_{out}$}
 * (3,0) to [C, l=$C$, *-] (3,-1.5)
 * node[ground]{$P_C$};
 * \draw [-latex] (5.2,0) -- (6,0) node[right] {$Q_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * # Background
 *
 * Zero-dimensional (0D) models
 * are lightweight methods to simulate bulk hemodynamic quantities in the
 * cardiovascular system. Unlike 3D and 1D models, 0D models are purely
 * time-dependent; they are unable to simulate spatial patterns in the
 * hemodynamics. 0D models are analogous to electrical circuits. The flow rate
 * simulated by 0D models represents electrical current, while the pressure
 * represents voltage. Three primary building blocks of 0D models are resistors,
 * capacitors, and inductors Resistance captures the viscous effects of blood
 * flow, capacitance represents the compliance and distensibility of the vessel
 * wall, and inductance represents the inertia of the blood flow. Different
 * combinations of these building blocks, as well as others, can be formed to
 * reflect the hemodynamics and physiology of different cardiovascular
 * anatomies.These 0D models are governed by differential algebraic equations
 * (DAEs).
 *
 * # Architecture
 *
 * ## Model
 *
 * The solver uses a highly modular framework to model the vascular anatomy,
 * using individual classes to represent different 0D elements of the
 * model. The elements are part of the MODEL namespace. Currently
 * supported elements are:
 *
 * * MODEL::BloodVessel: RCL blood vessel respresentation with optional
 * stenosis.
 * * MODEL::Junction: Junction element with arbitrary inlets and outlets.
 * * MODEL::FlowReferenceBC: Prescribed flow boundary condition.
 * * MODEL::PressureReferenceBC: Prescribed pressure boundary condition.
 * * MODEL::ResistanceBC: Resistance boundary condition.
 * * MODEL::WindkesselBC: RCR Windkessel boundary condition.
 * * MODEL::OpenLoopCoronaryBC: Open Loop coronary boundary condition.
 *
 * The elements are based on the parent MODEL::Block class. More information
 * about the elements can be found on their respective pages. The elements are
 * connected to each other via nodes (see MODEL::Node). Each node corresponds to
 * a flow and a pressure value of the model. The MODEL::DOFHandler handles the
 * degrees-of-freedom (DOF) of the system by assigning DOF indices to each
 * element that determine the location of the local element contribution in the
 * global system. The complete model is stored and managed by the MODEL::Model
 * class.
 *
 * ## Algebra
 *
 * Everything related to solving the system of equation is contained in the
 * ALGEBRA namespace. The svZeroDSolver provides two different strategies for
 * handling the system of equation: a dense strategy (see ALGEBRA::DenseSystem)
 * and a sparse strategy (ALGEBRA::SparseSystem). The ALGEBRA::Integrator
 * handles the time integration scheme (see Integrator documentation for
 * more information on that).
 *
 * ## Input/output
 *
 * Finally the IO namespace provides everything related to reading and writing
 * files. The IO::ConfigReader reads the simultion config and IO::write_csv and
 * IO::write_json are two different methods for writing the result files.
 *
 * All classes make use of templates to allow easy exchange of scalar
 * types like `double` or `float` or to switch between sparse and dense systems.
 *
 * # Build svZeroDSolver
 *
 * svZeroDSolver uses CMake to build the tool. The dependencies are installed
 * by CMake. If you want to use the C++ Python interface of svZeroDSolver make
 * sure to checkout the correct Python environment beforehand. Currently
 * supported are:
 * 
 * * Ubuntu 18
 * * Ubuntu 20
 * * Ubuntu 22
 * * macOS Big Sur
 * * macOS Monterey
 *
 * ## Build in debug mode
 *
 * \code
 * mkdir Debug
 * cd Debug
 * cmake -DCMAKE_BUILD_TYPE=Debug ..
 * cmake --build .
 * \endcode
 *
 * ## Build in release mode
 * \code
 * mkdir Release
 * cd Release
 * cmake -DCMAKE_BUILD_TYPE=Release ..
 * cmake --build .
 * \endcode
 * 
 * ## Build on Sherlock
 * \code
 * module load cmake/3.23.1 gcc/12.1.0 binutils/2.38
 * mkdir Release
 * cd Release
 * cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/share/software/user/open/gcc/12.1.0/bin/g++ -DCMAKE_C_COMPILER=/share/software/user/open/gcc/12.1.0/bin/gcc ..
 * cmake --build .
 * \endcode
 *
 * # Run svZeroDSolver
 *
 * After building svZeroDSolver, the build folder contains an executable
 * called `svzerodsolver`. Run it with one of:
 *
 * \code
 * ./svzerodsolver path/to/config.json path/to/output.csv
 * \endcode
 *
 * `path/to/config.json` and `path/to/output.csv` should be replaced by the
 * correct paths to the input and output file, respectively.
 *
 * ## Simulation parameters
 *
 * The svZeroDSolver can be configured with the following options in the
 * `simulation_parameters` section of the input file. Parameters without a
 * default value must be specified.
 * 
 *
 * Parameter key                           | Description                               | Default value
 * --------------------------------------- | ----------------------------------------- | ----------- 
 * number_of_cardiac_cycles                | Number of cardiac cycles to simulate      | -
 * number_of_time_pts_per_cardiac_cycle    | Number of time steps per cardiac cycle    | -
 * absolute_tolerance                      | Absolute tolerance for time integration   | \f$10^{-8}\f$
 * maximum_nonlinear_iterations            | Maximum number of nonlinear iterations for time integration | \f$30\f$
 * steady_initial                          | Toggle whether to use the steady solution as the initial condition for the simulation | true
 * output_variable_based                   | Output solution based on variables (i.e. flow+pressure at nodes and internal variables) | false
 * output_interval                         | The frequency of writing timesteps to the output (1 means every time step is written to output) | \f$1\f$
 * output_mean_only                        | Write only the mean values over every timestep to output file | false
 * output_derivative                       | Write time derivatives to output file | false
 * output_last_cycle_only                  | Write only the last cardiac cycle | false
 * 
 * 
 * # Developer Guide
 * 
 * If you are a developer and want to contribute to svZeroDSolver, you can find
 * more helpful information in our [Developer Guide](docs/cpp/pages/developer_guide.md).
 * 
 */
// clang-format on
#include "algebra/integrator.hpp"
#include "algebra/state.hpp"
#include "helpers/debug.hpp"
#include "helpers/endswith.hpp"
#include "io/configreader.hpp"
#include "io/csvwriter.hpp"
#include "model/model.hpp"

typedef double T;

template <typename TT>
using S = ALGEBRA::SparseSystem<TT>;

/**
 *
 * @brief Run svZeroDSolver with configuration
 *
 * 1. Read the input file
 * 2. Create the 0D model
 * 3. (Optional) Solve for steady initial condition
 * 4. Run simulation
 * 5. Write output to file
 *
 * @param json_config Path config or json encoded string with config
 * @return Result as csv encoded string
 */
const std::string run(std::string& json_config) {
  // Create configuration reader
  IO::ConfigReader<T> config(json_config);

  // Create model
  DEBUG_MSG("Creating model");
  auto model = config.get_model();
  DEBUG_MSG("Size of system:      " << model.dofhandler.size());

  // Get simulation parameters
  DEBUG_MSG("Setup simulutation");
  T time_step_size = config.get_time_step_size();
  DEBUG_MSG("Time step size:      " << time_step_size);
  int num_time_steps = config.get_num_time_steps();
  DEBUG_MSG("Number of timesteps: " << num_time_steps);
  T absolute_tolerance =
      config.get_scalar_simulation_parameter("absolute_tolerance", 1e-8);
  int max_nliter =
      config.get_int_simulation_parameter("maximum_nonlinear_iterations", 30);
  bool steady_initial =
      config.get_bool_simulation_parameter("steady_initial", true);
  bool output_variable_based =
      config.get_bool_simulation_parameter("output_variable_based", false);
  int output_interval =
      config.get_int_simulation_parameter("output_interval", 1);
  bool output_mean_only =
      config.get_bool_simulation_parameter("output_mean_only", false);
  bool output_derivative =
      config.get_bool_simulation_parameter("output_derivative", false);
  bool output_last_cycle_only =
      config.get_bool_simulation_parameter("output_last_cycle_only", false);

  // Setup system
  DEBUG_MSG("Starting simulation");
  ALGEBRA::State<T> state = ALGEBRA::State<T>::Zero(model.dofhandler.size());

  // Create steady initial
  if (steady_initial) {
    DEBUG_MSG("Calculating steady initial condition");
    T time_step_size_steady = config.cardiac_cycle_period / 10.0;
    auto model_steady = config.get_model();
    model_steady.to_steady();
    ALGEBRA::Integrator<T, S> integrator_steady(model_steady,
                                                time_step_size_steady, 0.1,
                                                absolute_tolerance, max_nliter);
    for (size_t i = 0; i < 31; i++) {
      state = integrator_steady.step(state, time_step_size_steady * T(i),
                                     model_steady);
    }
  }

  ALGEBRA::Integrator<T, S> integrator(model, time_step_size, 0.1,
                                       absolute_tolerance, max_nliter);

  std::vector<ALGEBRA::State<T>> states;
  std::vector<T> times;
  states.reserve(num_time_steps);
  times.reserve(num_time_steps);

  T time = 0.0;

  states.push_back(state);
  times.push_back(time);

  int interval_counter = 0;
  for (size_t i = 1; i < num_time_steps; i++) {
    state = integrator.step(state, time, model);
    interval_counter += 1;
    time = time_step_size * T(i);
    if (interval_counter == output_interval) {
      times.push_back(time);
      states.push_back(std::move(state));
      interval_counter = 0;
    }
  }
  DEBUG_MSG("Simulation completed");

  // Extract last cardiac cycle
  if (output_last_cycle_only) {
    states.erase(states.begin(), states.end() - config.num_pts_per_cycle);
    times.erase(times.begin(), times.end() - config.num_pts_per_cycle);
    T start_time = times[0];
    for (auto& time : times) {
      time -= start_time;
    }
  }

  std::string output;
  if (output_variable_based) {
    output = IO::to_variable_csv<T>(times, states, model, output_mean_only,
                                    output_derivative);
  } else {
    output = IO::to_vessel_csv<T>(times, states, model, output_mean_only,
                                  output_derivative);
  }
  return output;
}
