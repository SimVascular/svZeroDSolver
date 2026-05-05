@page calibrator svZeroDCalibrator

[TOC]

svZeroDCalibrator solves the inverse problem: given observations of pressure
and flow, it tunes the parameters of an existing 0D model so that the
solver-predicted state matches the observations. The optimizer is a
Levenberg-Marquardt loop running over the alpha vector of all calibratable
parameters in the model. See `src/optimize/calibrate.cpp` and
`src/optimize/LevenbergMarquardtOptimizer.h` for the implementation.

## Input file

The input file is a normal solver configuration extended with three extra
blocks:

```json
{
  "boundary_conditions": [...],
  "junctions": [...],
  "vessels": [...],
  "calibration_parameters": { ... },
  "y": { ... },
  "dy": { ... }
}
```

* `y` and `dy` are dictionaries keyed by the variable names that the model
  exposes after assembly (e.g. `pressure:INFLOW:branch0_seg0`,
  `flow:branch0_seg0:OUT`). Each value is a list of observed samples; the
  number of samples is the same for every variable. `dy` holds the matching
  time derivatives.
* `calibration_parameters` collects the calibrator-specific options described
  below.

The output file has the same shape as a solver input but with calibrated
values written into `zero_d_element_values` (vessels) and `junction_values`
(multi-outlet junctions). The `y`, `dy` and `calibration_parameters` keys are
removed from the output.

## Selecting which parameters to calibrate

By default every parameter exposed by every supported block is calibrated.
Two optional fields restrict the set of parameters that are optimized; any
parameter that is not selected is held constant at the value found in the
input file.

### Global default: `calibration_parameters.calibrate`

The list of parameter names under `calibration_parameters.calibrate` applies
to every block that does not override it. The names must match the parameter
names a block exposes through its `input_params` list (e.g. `R_poiseuille`,
`C`, `L`, `stenosis_coefficient` for `BloodVessel`).

```json
"calibration_parameters": {
  "tolerance_gradient": 1e-5,
  "tolerance_increment": 1e-10,
  "maximum_iterations": 100,
  "initial_damping_factor": 1.0,
  "calibrate": ["R_poiseuille"]
}
```

In the example above only the Poiseuille resistance is calibrated for every
block; capacitance, inductance, and stenosis coefficient stay at the values
provided by the input file.

### Per-block override: `calibrate` field on a vessel or junction

A vessel or junction can carry its own `calibrate` field at the top level of
its block entry. When present, this list takes precedence over the global
default for that one block:

```json
{
  "vessel_id": 0,
  "vessel_name": "branch0_seg0",
  "zero_d_element_type": "BloodVessel",
  "zero_d_element_values": {
    "R_poiseuille": 0.0,
    "C": 1.2e-6,
    "L": 0.25,
    "stenosis_coefficient": 1.06e-5
  },
  "calibrate": ["R_poiseuille"]
}
```

```json
{
  "junction_name": "J0",
  "junction_type": "BloodVesselJunction",
  "junction_values": {
    "R_poiseuille": [0.0, 0.0],
    "L": [0.0, 0.0],
    "stenosis_coefficient": [0.0, 0.0]
  },
  "inlet_vessels": [0],
  "outlet_vessels": [1, 2],
  "calibrate": ["R_poiseuille"]
}
```

### Resolution order

For each block, the calibrator picks the active set of parameter names with
the following precedence:

1. The block's own `calibrate` field, if present.
2. Otherwise, `calibration_parameters.calibrate`.
3. Otherwise, every parameter the block exposes (legacy behavior).

An explicit empty list (`"calibrate": []`) at any level means "calibrate
nothing for that scope". The calibrator errors out if no parameter is
selected in any block.

### Worked examples

The repository ships three small fixtures derived from the
`0104_0001` Vascular Model Repository case that exercise each path:

* `tests/cases/vmr/input/0104_0001_calibrate_R_only_global.json` -
  uses `calibration_parameters.calibrate`.
* `tests/cases/vmr/input/0104_0001_calibrate_R_only_per_block.json` -
  uses per-block `calibrate` fields with no global default.
* `tests/cases/vmr/input/0104_0001_calibrate_R_only_block_overrides.json` -
  sets a misleading global default and overrides every block; the override
  must win.

In each case every `R_poiseuille` value is zeroed in the input file while
`C`, `L`, and `stenosis_coefficient` are kept at their ground-truth values.
The calibrator recovers the reference R values to machine precision and
leaves the other parameters untouched. The matching test is
`tests/test_calibrator.py::test_calibration_R_only`.

## Block requirements

A block is calibratable as long as it implements `update_gradient` and
exposes its parameter names through the standard `Block::input_params`
field. The calibrator reads this metadata at runtime via
`Block::input_params` and `Block::input_params_list`, so adding a new
calibratable block does not require any changes to `calibrate.cpp`.

The legacy flag `calibration_parameters.calibrate_stenosis_coefficient`
(default `true`) layers on top of the selection logic: when set to `false`,
`stenosis_coefficient` is held constant regardless of any `calibrate` field.
The flag predates the `calibrate` field and is preserved for backward
compatibility; new input files should prefer `calibrate`.
