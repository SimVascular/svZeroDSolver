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
 * * <a href="https://github.com/SimVascular/svZeroDSolver">Source repository</a>
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
 * cardiovascular system. Unlike 3D and 1D models, 0D models are purely time-dependent;
 * they are unable to simulate spatial patterns in the hemodynamics.
 * 0D models are analogous to electrical circuits. The flow rate simulated by
 * 0D models represents electrical current, while the pressure represents voltage.
 * Three primary building blocks of 0D models are resistors, capacitors, and inductors
 * Resistance captures the viscous effects of blood flow, capacitance represents
 * the compliance and distensibility of the vessel wall, and inductance represents
 * the inertia of the blood flow. Different combinations of these building blocks,
 * as well as others, can be formed to reflect the hemodynamics and physiology of
 * different cardiovascular anatomies.These 0D models are governed
 * by differential algebraic equations (DAEs).
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
 * * MODEL::BloodVessel: RCL blood vessel respresentation with optional stenosis
 * * MODEL::Junction: Junction element with arbitrary inlets and outlets.
 * * MODEL::FlowReference: Prescribed flow boundary condition
 * * MODEL::RCRBlockWithDistalPressure: RCR Windkessel boundary condition
 *
 * More information about the elements can be found on their
 * respective pages. The elements are connected to each other via nodes
 * (see MODEL::Node). Each node corresponds to a flow and a pressure value
 * of the model. The MODEL::DOFHandler handles the degrees-of-freedom (DOF) of the
 * system by assigning DOF indices to each element that determine the location
 * of the local element contribution in the global system. The complete model
 * is stored and managed by the MODEL::Model class.
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
 * svZeroDSolver can be build easily via CMake. Make sure you have the following
 * dependencies installed before you start:
 *
 * ## Install dependencies on macOS
 *
 * \code
 * brew install eigen       # Linear algebra library
 * brew install jsoncpp     # Standard json library
 * brew install simdjson    # Fast json input parser
 * \endcode
 *
 * ## Install dependencies on Linux
 *
 * \code
 * sudo apt install libeigen3-dev       # Linear algebra library
 * sudo apt-get install libjsoncpp-dev  # Standard json library
 * sudo apt install libsimdjson-dev     # Fast json input parser
 * \endcode
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
 * # Run svZeroDSolver
 *
 * After building svZeroDSolver, the build folder contains an executable
 * called `svzerodsolver`. Run it with one of:
 *
 * \code
 * ./svzerodsolver path/to/config.json path/to/output.json  # For json output file format
 * ./svzerodsolver path/to/config.json path/to/output.csv   # For csv output file format (faster)
 * \endcode
 *
 * `path/to/config.json` and `path/to/output` should be replaced by the correct
 * paths to the input and output file, respectively.
 *
 */