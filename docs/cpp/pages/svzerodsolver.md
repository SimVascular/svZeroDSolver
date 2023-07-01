@page svzerodsolver svZeroDSolver

[TOC]

# Run svZeroDSolver

## From the command line

svZeroDSolver can be executed from the command line using a JSON configuration
file.

```bash
svzerodsolver tests/cases/steadyFlow_RLC_R.json result_steadyFlow_RLC_R.csv
```

The result will be written to a csv file.

## In Python

For some applications it is beneficial to run svZeroDSolver directly
from within another Python application. This can, for example, be
useful when many simulations need to be performed (e.g. for 
calibration, uncertainty quantification, ...). Please make sure that
you installed svZerodPlus via pip to enable this feature. We start by
importing svzerodplus:

```python
import svzerodplus
```

Next, we create a solver from our configuration. The configuration can
be specified by either a path to a JSON file:

```python
solver = svzerodplus.Solver("tests/cases/steadyFlow_RLC_R.json")
```

or as a Python dictionary:

```python
my_config = {...}
solver = svzerodplus.Solver(my_config)
```

To run the simulation we add:

```python
solver.run()
```

The simulation result is now saved in the solver instance. We can obtain
results for individual degrees-of-freedom (DOFs) as
```python
solver.get_single_result("flow:INFLOW:branch0_seg0")

>>> array([5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.,
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
solver.get_single_result_avg("flow:INFLOW:branch0_seg0")

>>> 5.0
```


# Configuration (file)

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

## Simulation parameters

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


## Vessels

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


## Junctions

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

## Boundary conditions

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
