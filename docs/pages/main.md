@mainpage svZeroDPlus

[TOC]

svZeroDPlus is an application for performing simulations with 0D/lumped-parameter
computer models for cardiovascular flows. 

Some noteworthy features of svZeroDPlus are:
* It is completely modular. Users can create custom flow models by arranging 
blocks corresponding to blood vessels, junctions, different types of 
boundary conditions, etc. 
* It is written in C++ to enable high-performance applications.
* svZeroDPlus offers both a Python API and a C++ shared library to interface with other 
Python or C++-based applications. This allows it to be used in a fully coupled manner 
with other multi-physics solvers, and for parameter estimation, uncertainty 
quantification, etc.
* The svZeroDCalibrator application, which is included in svZeroDPlus, optimizes 0D
blood vessel parameters to recapitulate given time-varying flow and pressure measurements 
(for example, from a high-fidelity 3D simulation). This allows users to build accurate 
0D models that reflect observed hemodynamics.

Zero-dimensional (0D) models
are lightweight methods to simulate bulk hemodynamic quantities in the
cardiovascular system. Unlike 3D and 1D models, 0D models are purely
time-dependent; they are unable to simulate spatial patterns in the
hemodynamics. 0D models are analogous to electrical circuits. The flow rate
simulated by 0D models represents electrical current, while the pressure
represents voltage. Three primary building blocks of 0D models are resistors,
capacitors, and inductors Resistance captures the viscous effects of blood
flow, capacitance represents the compliance and distensibility of the vessel
wall, and inductance represents the inertia of the blood flow. Different
combinations of these building blocks, as well as others, can be formed to
reflect the hemodynamics and physiology of different cardiovascular
anatomies. These 0D models are governed by differential algebraic equations
(DAEs).

For more background information on 0D models, have a look at SimVascular's
[ROM Simulation Guide](http://simvascular.github.io/docsROMSimulation.html).

* <a href="https://github.com/StanfordCBCL/svZeroDPlus">Source
repository</a>
* <a href="https://simvascular.github.io">About SimVascular</a>

# Installation

svZeroDPlus can be installed in two different ways. For using the Python
API, an installation via pip is recommended.

## Using pip

For a pip installation, simply run the following command
(cloning of the repository is not required):

```bash
pip install git+https://github.com/StanfordCBCL/svZeroDPlus.git
```

## Using CMake

If you want to build svZeroDPlus manually from source, clone the repository
and run the following commands from the top directory of the project:

```bash
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```


@remark <details>
  <summary>**Building on Sherlock**</summary>

```bash
module load cmake/3.23.1 gcc/12.1.0 binutils/2.38
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/share/software/user/open/gcc/12.1.0/bin/g++ -DCMAKE_C_COMPILER=/share/software/user/open/gcc/12.1.0/bin/gcc ..
cmake --build .
```

</details>

# Developer Guide

If you are a developer and want to contribute to svZeroDSolver, you can find
more helpful information in our [Developer Guide](@ref developer_guide).

# svZeroDSolver

svZeroDSolver can be used to run zero-dimensional (0D) cardiovascular
simulations based on a given configuration.

## Run svZeroDSolver from the command line

svZeroDSolver can be executed from the command line using a JSON configuration
file.

```bash
svzerodsolver tests/cases/steadyFlow_RLC_R.json result_steadyFlow_RLC_R.csv
```

The result will be written to a csv file.

## Run svZeroDSolver from other programs

For some applications it is beneficial to run svZeroDSolver directly
from within another program. For example, this can be
useful when many simulations need to be performed (e.g. for 
calibration, uncertainty quantification, ...). It is also allows using
svZeroDPlus with other solvers, for example as boundary conditions or
forcing terms.

### In C++

SvZeroDPlus needs to be built using CMake to use the shared library interface.

Detailed examples of interfacing with svZeroDPlus from C++ codes are available 
in the test cases at `svZeroDPlus/tests/test_interface`. 

### In Python

Please make sure that
you installed svZerodPlus via pip to enable this feature. We start by
importing svzerodplus:

```python
>>> import svzerodplus
```

Next, we create a solver from our configuration. The configuration can
be specified by either a path to a JSON file:

```python
>>> solver = svzerodplus.Solver("tests/cases/steadyFlow_RLC_R.json")
```

or as a Python dictionary:

```python
>>> my_config = {...}
>>> solver = svzerodplus.Solver(my_config)
```

To run the simulation we add:

```python
>>> solver.run()
```

The simulation result is now saved in the solver instance. We can obtain
results for individual degrees-of-freedom (DOFs) as
```python
>>> solver.get_single_result("flow:INFLOW:branch0_seg0")

array([5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.,
       5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.,
       5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.,
       5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.,
       5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.,
       5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.])
```
The naming of the DOFs is similar to how results are written if the simulation
option `output_variable_based` is activated (see below). We can also obtain
the mean result for a DOF over time with:
```python
>>> solver.get_single_result_avg("flow:INFLOW:branch0_seg0")

5.0
```

Or the result of the full simulation as a pandas data frame:

```python
>>> solver.get_full_result()

             name  time  flow_in  flow_out  pressure_in  pressure_out
0    branch0_seg0  0.00      5.0       5.0       1100.0         600.0
1    branch0_seg0  0.01      5.0       5.0       1100.0         600.0
2    branch0_seg0  0.02      5.0       5.0       1100.0         600.0
3    branch0_seg0  0.03      5.0       5.0       1100.0         600.0
4    branch0_seg0  0.04      5.0       5.0       1100.0         600.0
..            ...   ...      ...       ...          ...           ...
96   branch0_seg0  0.96      5.0       5.0       1100.0         600.0
97   branch0_seg0  0.97      5.0       5.0       1100.0         600.0
98   branch0_seg0  0.98      5.0       5.0       1100.0         600.0
99   branch0_seg0  0.99      5.0       5.0       1100.0         600.0
100  branch0_seg0  1.00      5.0       5.0       1100.0         600.0

[101 rows x 6 columns]
```

There is also a function to retrieve the full result directly based on a given configuration:

```python

>>> my_config = {...}
>>> svzerodplus.simulate(my_config)

             name  time  flow_in  flow_out  pressure_in  pressure_out
0    branch0_seg0  0.00      5.0       5.0       1100.0         600.0
1    branch0_seg0  0.01      5.0       5.0       1100.0         600.0
2    branch0_seg0  0.02      5.0       5.0       1100.0         600.0
3    branch0_seg0  0.03      5.0       5.0       1100.0         600.0
4    branch0_seg0  0.04      5.0       5.0       1100.0         600.0
..            ...   ...      ...       ...          ...           ...
96   branch0_seg0  0.96      5.0       5.0       1100.0         600.0
97   branch0_seg0  0.97      5.0       5.0       1100.0         600.0
98   branch0_seg0  0.98      5.0       5.0       1100.0         600.0
99   branch0_seg0  0.99      5.0       5.0       1100.0         600.0
100  branch0_seg0  1.00      5.0       5.0       1100.0         600.0

[101 rows x 6 columns]

```


## Configuration

svZeroDSolver is configured using either a JSON file or a Python
dictionary. The top-level structure of both is:

```python
{
    "simulation_parameters": {...},
    "vessels": [...],
    "junctions": [...],
    "boundary_conditions": [...]
}
```

In the following sections, the individual categories are described in more
detail.

### Simulation parameters

The svZeroDSolver can be configured with the following options in the
`simulation_parameters` section of the input file. Parameters without a
default value must be specified.


Parameter key                           | Description                               | Default value
--------------------------------------- | ----------------------------------------- | ----------- 
number_of_cardiac_cycles                | Number of cardiac cycles to simulate      | -
number_of_time_pts_per_cardiac_cycle    | Number of time steps per cardiac cycle    | -
absolute_tolerance                      | Absolute tolerance for time integration   | \f$10^{-8}\f$
maximum_nonlinear_iterations            | Maximum number of nonlinear iterations for time integration | \f$30\f$
steady_initial                          | Toggle whether to use the steady solution as the initial condition for the simulation | true
output_variable_based                   | Output solution based on variables (i.e. flow+pressure at nodes and internal variables) | false
output_interval                         | The frequency of writing timesteps to the output (1 means every time step is written to output) | \f$1\f$
output_mean_only                        | Write only the mean values over every timestep to output file | false
output_derivative                       | Write time derivatives to output file | false
output_all_cycles                       | Write all cardiac cycles to output file | false


### Vessels

More information about the vessels can be found in their respective class references.

```python
{
    "boundary_conditions": {
        "inlet": "INFLOW", # Optional: Name of inlet boundary condition
        "outlet": "OUT", # Optional: Name of outlet boundary condition
    },
    "vessel_id": 0, # ID of the vessel
    "vessel_name": "branch0_seg0", # Name of vessel
    "zero_d_element_type": "BloodVessel", # Type of vessel
    "zero_d_element_values": {...} # Values for configuration parameters
}
```

Description                           | Class                       | `zero_d_element_type` | `zero_d_element_values`
------------------------------------- | --------------------------- | --------------------- | ------------------------
Blood vessel with \n optional stenosis   | MODEL::BloodVessel       | `BloodVessel`         | `C`: Capacity \n `L`: Inductance \n `R_poiseuille`: Poiseuille resistance \n `stenosis_coefficient`: Stenosis coefficient


### Junctions

More information about the junctions can be found in their respective class references.

```python
{
    "junction_name": "J0", # Name of the junction
    "junction_type": "BloodVesselJunction", # Type of the junction
    "inlet_vessels": [0], # List of vessel IDs connected to the inlet
    "outlet_vessels": [1, 2], # List of vessel IDs connected to the inlet
    "junction_values": {...} # Values for configuration parameters
}
```

Description                           | Class                       | `junction_type`       | `junction_values`
------------------------------------- | --------------------------- | --------------------- | ----------- 
Purely mass \n conserving \n junction | MODEL::Junction             | `NORMAL_JUNCTION`     | -
Resistive \n junction                 | MODEL::ResistiveJunction    | `resistive_junction`  | `R`: Ordered list of resistances for all inlets and outlets
Blood vessel \n junction              | MODEL::BloodVesselJunction  | `BloodVesselJunction` | Same as for `BloodVessel` element but \n as ordered list for each inlet and outlet

### Boundary conditions

```python
{
    "bc_name": "INFLOW", # Name of the boundary condition
    "bc_type": "FLOW", # Type of the boundary condition
    "bc_values": {...} # Values for configuration parameters
},
```

Description                           | Class                       | `bc_type`             | `bc_values`
------------------------------------- | --------------------------- | --------------------- | ----------- 
Prescribed (transient) flow           | MODEL::FlowReferenceBC      | `FLOW`                | `Q`: Time-dependent flow values \n `t`: Time stamps
Prescribed (transient) pressure       | MODEL::PressureReferenceBC  | `PRESSURE`            | `P`: Time-dependent pressure values \n `t`: Time stamps
Resistance                            | MODEL::ResistanceBC         | `RESISTANCE`          | `R`: Resistance \n `Pd`: Time-dependent distal pressure \n `t`: Time stamps
Windkessel                            | MODEL::WindkesselBC         | `RCR`                 | `Rp`: Proximal resistance \n `C`: Capacitance \n `Rd`: Distal resistance \n `Pd`: Distal pressure


# svZeroDCalibrator

svZeroDCalibrator can be used to calibrate cardiovascular 0D models (i.e. infer optimal
parameters for the 0D elements) based on a given transient result (i.e. from a
3D simulation).

## Run svZeroDCalibrator

### From the command line
svZeroDCalibrator can be executed from the command line using a JSON configuration
file.

```bash
svzerodcalibrator path/to/input_file.json path/to/output_file.json
```

The result will be written to a JSON file.


### In Python

svZeroDCalibrator can also be called directly from Python.
Please make sure that you installed svZerodPlus via pip to enable this feature. We start by
importing svzerodplus:

```python
import svzerodplus

my_unoptimized_config = {...}
my_optimized_config = svzerodplus.calibrate(my_unoptimized_config)
```

## Configuration (file)

In order to make svZeroDCalibrator easy to use, it is based on a similar configuration
file than svZeroDSolver. Instead of the `simulation_parameters` section, it has a section
called `calibration_parameters`. Additionally the optimization target (i.e. a given)
3D result is passed with the key `y` and it's temporal derivative via `dy`. See
`tests/cases/steadyFlow_calibration.json` for an example input file.

```python
{
    "calibration_parameters": {...},
    "vessels": [...],
    "junctions": [...],
    "boundary_conditions": [...],
    "y": {
      "flow:INFLOW:branch0_seg0": [0.0, 0.1, ...],  # Time series for DOF
      "pressure:INFLOW:branch0_seg0": [0.0, 0.1, ...],  # Time series for DOF
      ...
    },
    "dy": {
      "flow:INFLOW:branch0_seg0": [0.0, 0.1, ...],  # Time series for DOF
      "pressure:INFLOW:branch0_seg0": [0.0, 0.1, ...],  # Time series for DOF
      ...
    },
}
```

### Calibration parameters

Here is a list of the parameters that can be specified in the `calibration_parameters`
section of the configuration file.

Parameter key                           | Description                               | Default value
--------------------------------------- | ----------------------------------------- | ----------- 
tolerance_gradient                      | Gradient tolerance for calibration        | \f$10^{-5}\f$
tolerance_increment                     | Increment tolerance for calibration       | \f$10^{-10}\f$
maximum_iterations                      | Maximum calibration iterations            | 100
calibrate_stenosis_coefficient          | Toggle whether stenosis coefficient should be calibrated        | True
set_capacitance_to_zero                 | Toggle whether all capacitances should be manually set to zero  | False
initial_damping_factor                  | Initial damping factor for Levenberg-Marquardt optimization  | 1.0
