@mainpage svZeroDSolver

[TOC]

svZeroDSolver is an application for performing simulations with 0D/lumped-parameter
computer models for cardiovascular flows. 

Some noteworthy features of svZeroDSolver are:
* It is completely modular. Users can create custom flow models by arranging 
blocks corresponding to blood vessels, junctions, different types of 
boundary conditions, etc. 
* It is written in C++ to enable high-performance applications.
* svZeroDSolver offers both a Python API and a C++ shared library to interface with other 
Python or C++-based applications. This allows it to be used in a fully coupled manner 
with other multi-physics solvers, and for parameter estimation, uncertainty 
quantification, etc.
* The svZeroDCalibrator application, which is included in svZeroDSolver, optimizes 0D
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
capacitors, and inductors. Resistance captures the viscous effects of blood
flow, capacitance represents the compliance and distensibility of the vessel
wall, and inductance represents the inertia of the blood flow. Different
combinations of these building blocks, as well as others, can be formed to
reflect the hemodynamics and physiology of different cardiovascular
anatomies. These 0D models are governed by differential algebraic equations
(DAEs).

The main categories of blocks implemented in svZeroDSolver are:
- Blood vessels
- Junctions
- Boundary conditions
- Heart chambers
- Heart valves

For an overview of available 0D elements (blocks) see: Block

You can find more details about governing equations in individual blocks in their respective documentation pages. For example:
- BloodVessel
- BloodVesselJunction
- WindkesselBC

For implementation details, have a look at the [source code](https://github.com/simvascular/svZeroDSolver).
Mathematics details can be found in the following classes:
- System of equations: SparseSystem
- Time integration: Integrator

[More information about SimVascular](https://simvascular.github.io)

# Installation

svZeroDSolver can be installed in two different ways. For using the Python
API, an installation via pip is recommended.

## Using pip

For a pip installation, simply run the following command
(cloning of the repository is not required):

```bash
pip install git+https://github.com/simvascular/svZeroDSolver.git
```

## Using CMake

If you want to build svZeroDSolver manually from source, clone the repository
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

# svZeroDSolver - Quick User Guide

svZeroDSolver can be used to run zero-dimensional (0D) cardiovascular
simulations based on a given configuration.

## Run svZeroDSolver from the command line

svZeroDSolver can be executed from the command line using a JSON configuration
file.

```bash
svzerodsolver tests/cases/steadyFlow_RLC_R.json result_steadyFlow_RLC_R.csv
```

The result will be written to a CSV file.

## Run svZeroDSolver from other programs

For some applications it is beneficial to run svZeroDSolver directly
from within another program. For example, this can be
useful when many simulations need to be performed (e.g. for 
calibration, uncertainty quantification, ...). It is also allows using
svZeroDSolver with other solvers, for example as boundary conditions or
forcing terms.

### In C++

SvZeroDSolver needs to be built using CMake to use the shared library interface.

Detailed examples of interfacing with svZeroDSolver from C++ codes are available 
in the test cases at `svZeroDSolver/tests/test_interface`. 

### In Python

Please make sure that
you installed svZerodSolver via pip to enable this feature. We start by
importing pysvzerod:

```python
>>> import pysvzerod
```

Next, we create a solver from our configuration. The configuration can
be specified by either a path to a JSON file:

```python
>>> solver = pysvzerod.Solver("tests/cases/steadyFlow_RLC_R.json")
```

or as a Python dictionary:

```python
>>> my_config = {...}
>>> solver = pysvzerod.Solver(my_config)
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
>>> pysvzerod.simulate(my_config)

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


Parameter key                             | Description                               | Default value
----------------------------------------- | ----------------------------------------- | ----------- 
`number_of_cardiac_cycles`                | Number of cardiac cycles to simulate      | -
`number_of_time_pts_per_cardiac_cycle`    | Number of time steps per cardiac cycle    | -
`absolute_tolerance`                      | Absolute tolerance for time integration   | \f$10^{-8}\f$
`maximum_nonlinear_iterations`            | Maximum number of nonlinear iterations for time integration | \f$30\f$
`steady_initial`                          | Toggle whether to use the steady solution as the initial condition for the simulation | true
`output_variable_based`                   | Output solution based on variables (i.e. flow+pressure at nodes and internal variables) | false
`output_interval`                         | The frequency of writing timesteps to the output (1 means every time step is written to output) | \f$1\f$
`output_mean_only`                        | Write only the mean values over every timestep to output file | false
`output_derivative`                       | Write time derivatives to output file | false
`output_all_cycles`                       | Write all cardiac cycles to output file | false
`use_cycle_to_cycle_error`                | Use cycle-to-cycle error to determine number of cycles for convergence | false
`sim_cycle_to_cycle_percent_error`        | Percentage error threshold for cycle-to-cycle pressure and flow difference | 1.0

The option `use_cycle_to_cycle_error` allows the solver to change the number of cardiac cycles it runs depending on the cycle-to-cycle convergence of the simulation. For simulations with no RCR boundary conditions, the simulation will add extra cardiac cycles until the difference between the mean pressure and flow in consecutive cycles is below the threshold set by `sim_cycle_to_cycle_percent_error` at all inlets and outlets of the model. If there is at least one RCR boundary condition, the number of cycles is determined based on equation 21 of \cite pfaller21, using the RCR boundary condition with the largest time constant.

### Vessels

More information about the vessels can be found in their respective class references. Below is a template vessel block with boundary conditions, `INFLOW` and `OUT`, at its inlet and outlet respectively.

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

Description                              | Class                       | `zero_d_element_type` | `zero_d_element_values`
---------------------------------------- | --------------------------- | --------------------- | ------------------------
Blood vessel with \n optional stenosis   | BloodVessel                 | `BloodVessel`         | `C`: Capacity \n `L`: Inductance \n `R_poiseuille`: Poiseuille resistance \n `stenosis_coefficient`: Stenosis coefficient


### Junctions

More information about the junctions can be found in their respective class references. Below is a template junction block that connects vessel ID 0 with vessel IDs 1 and 2.

```python
{
    "junction_name": "J0", # Name of the junction
    "junction_type": "BloodVesselJunction", # Type of the junction
    "inlet_vessels": [0], # List of vessel IDs connected to the inlet
    "outlet_vessels": [1, 2], # List of vessel IDs connected to the inlet
    "junction_values": {...} # Values for configuration parameters
}
```

Description                           | Class                | `junction_type`       | `junction_values`
------------------------------------- | ---------------------| --------------------- | ----------- 
Purely mass \n conserving \n junction | Junction             | `NORMAL_JUNCTION`     | -
Resistive \n junction                 | ResistiveJunction    | `resistive_junction`  | `R`: Ordered list of resistances for all inlets and outlets
Blood vessel \n junction              | BloodVesselJunction  | `BloodVesselJunction` | Same as for `BloodVessel` element but \n as ordered list for each inlet and outlet

### Boundary conditions

More information about the boundary conditions can be found in their respective class references. Below is a template `FLOW` boundary condition.

```python
{
    "bc_name": "INFLOW", # Name of the boundary condition
    "bc_type": "FLOW", # Type of the boundary condition
    "bc_values": {...} # Values for configuration parameters
},
```

Description                           | Class                  | `bc_type`             | `bc_values`
------------------------------------- | ---------------------- | --------------------- | ----------- 
Prescribed (transient) flow           | FlowReferenceBC        | `FLOW`                | `Q`: Time-dependent flow values \n `t`: Time steps \n `fn`: Mathematical expression \n Note: Either specify `Q` and `t` together, or just `fn`
Prescribed (transient) pressure       | PressureReferenceBC    | `PRESSURE`            | `P`: Time-dependent pressure values \n `t`: Time steps \n `fn`: Mathematical expression \n Note: Either specify `Q` and `t` together, or just `fn`
Resistance                            | ResistanceBC           | `RESISTANCE`          | `R`: Resistance \n `Pd`: Time-dependent distal pressure \n `t`: Time stamps
Windkessel                            | WindkesselBC           | `RCR`                 | `Rp`: Proximal resistance \n `C`: Capacitance \n `Rd`: Distal resistance \n `Pd`: Distal pressure
Coronary outlet                       | OpenLoopCoronaryBC     | `CORONARY`            | `Ra`: Proximal resistance \n `Ram`: Microvascular resistance \n `Rv`: Venous resistance \n `Ca`: Small artery capacitance \n `Cim`: Intramyocardial capacitance \n `Pim`: Intramyocardial pressure \n `Pv`: Venous pressure

The above table describes the most commonly used boundary conditions. In addition, svZeroDSolver includes various closed-loop boundary conditions. Examples can be found in `svZeroDSolver/tests/cases`.

Note that the `FLOW` and `PRESSURE` boundary conditions accept mathematical expressions in `bc_values`. For example, values of the boundary condition can be specified as a function of time as follow: 
```python
{
    "bc_name": "INFLOW", # Name of the boundary condition
    "bc_type": "FLOW", # Type of the boundary condition
    "bc_values": {
        "Q": [ ..., ..., ... ], # Comma-separated list of values
        "t": [ ..., ..., ... ]  # Comma-separated list of corresponding time stamps
    }
},
```
See `svZeroDSolver/tests/cases/pulsatileFlow_R_RCR.json` for an example.

They can also be specified as a mathematica expression as follow: 
```python
{
    "bc_name": "INFLOW", # Name of the boundary condition
    "bc_type": "FLOW", # Type of the boundary condition
    "bc_values": {
        "fn": "2.0 * (4*atan(1.)) * cos(2.0 * (4*atan(1.)) * t)"
    }
},
```
For an example with a mathematical expression for the boundary condition, see `svZeroDSolver/tests/cases/timeDep_Flow.json`. 

## Simulation Outputs

The siumulation outputs will be saved in the specified CSV file (`<name_of_output_file>.csv`) when running `svZeroDSolver` from the command line as follows:
```bash
svzerodsolver <name_of_configuration_file>.json <name_of_output_file>.csv
```
If the name of the CSV file is not specified, the default is `output.csv`. The format of the file depends on the user-specified configuration within the `simulation_parameters` block of the JSON configuration file. 

If `output_variable_based` is set to `true`, the CSV file will contain all the degrees-of-freedom in the simulation. Otherwise, only the flow and pressure at the inlets and outlets of vessels is written. 

The degrees-of-freedom (DOFs) follow the following naming scheme:

- Flow DOFs are labelled `flow:<name_of_upstream_block>:<name_of_downstream_block>`.
- Pressure DOFs are labelled `pressure:<name_of_upstream_block>:<name_of_downstream_block>`.
- Internal DOFs (i.e., variables internal to a block and not connected to upstream/downstream blocks) are labelled `<variable_name>:<block_name>`. The internal variables for each block are listed in the blocks' [class documentation](https://simvascular.github.io/svZeroDSolver/annotated.html). 

When the outputs are written in the variable-based and vessel-based forms, the user can specify whether they want outputs written for all cardiac cycles or just the last cardiac cycle using the `output_all_cycles` option. By default, only the last cycle is written. This makes the simulation more efficient. 

The number of timesteps between each time the output is written is specified by `output_interval`. By default, output is written at every time step. 


# svZeroDCalibrator - Quick User Guide

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
Please make sure that you installed svZerodSolver via pip to enable this feature. We start by
importing pysvzerod:

```python
import pysvzerod

my_unoptimized_config = {...}
my_optimized_config = pysvzerod.calibrate(my_unoptimized_config)
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
