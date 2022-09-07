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
 * to [R, l=$R$, *-] (3,0) node[anchor=south] {$P_{C}$}
 * (3,0) to [L, l=$L$, *-*] (5,0)
 * node[anchor=south]{$P_{out}$}
 * (3,0) to [C, l=$C$, *-] (3,-1.5)
 * node[ground]{};
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
 * ALGEBRA namespace. The ALGEBRA::Integrator
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
 * types like `double` or `float`.
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

// Setting scalar type to double
typedef double T;

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
  // Load model and configuration
  IO::ConfigReader<T> reader;
  reader.load(json_config);

  // Setup system
  ALGEBRA::State<T> state = reader.initial_state;

  // Create steady initial
  if (reader.sim_steady_initial) {
    T time_step_size_steady = reader.sim_cardiac_cycle_period / 10.0;
    reader.model.to_steady();
    ALGEBRA::Integrator<T> integrator_steady(
        reader.model, time_step_size_steady, 0.1, reader.sim_abs_tol,
        reader.sim_nliter);
    for (size_t i = 0; i < 31; i++) {
      state = integrator_steady.step(state, time_step_size_steady * T(i),
                                     reader.model);
    }
    reader.model.to_unsteady();
  }

  // Set-up integrator
  ALGEBRA::Integrator<T> integrator(reader.model, reader.sim_time_step_size,
                                    0.1, reader.sim_abs_tol, reader.sim_nliter);

  // Initialize loop
  std::vector<ALGEBRA::State<T>> states;
  std::vector<T> times;
  states.reserve(reader.sim_num_time_steps);
  times.reserve(reader.sim_num_time_steps);
  T time = 0.0;
  states.push_back(state);
  times.push_back(time);

  // Run integrator
  int interval_counter = 0;
  for (size_t i = 1; i < reader.sim_num_time_steps; i++) {
    state = integrator.step(state, time, reader.model);
    interval_counter += 1;
    time = reader.sim_time_step_size * T(i);
    if (interval_counter == reader.output_interval) {
      times.push_back(time);
      states.push_back(std::move(state));
      interval_counter = 0;
    }
  }

  // Extract last cardiac cycle
  if (reader.output_last_cycle_only) {
    states.erase(states.begin(), states.end() - reader.sim_pts_per_cycle);
    times.erase(times.begin(), times.end() - reader.sim_pts_per_cycle);
    T start_time = times[0];
    for (auto& time : times) {
      time -= start_time;
    }
  }

  // Write csv output string
  std::string output;
  if (reader.output_variable_based) {
    output = IO::to_variable_csv<T>(times, states, reader.model,
                                    reader.output_mean_only,
                                    reader.output_derivative);
  } else {
    output =
        IO::to_vessel_csv<T>(times, states, reader.model,
                             reader.output_mean_only, reader.output_derivative);
  }
  return output;
}
