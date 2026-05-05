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
* `calibration_parameters` collects the calibrator-specific options
  (tolerances, iteration cap, damping factor, capacitance clamp).

The output file has the same shape as a solver input but with calibrated
values written into `zero_d_element_values` (vessels) and `junction_values`
(multi-outlet junctions). The `y`, `dy` and `calibration_parameters` keys are
removed from the output.

## Selecting which parameters to calibrate

Every vessel and multi-outlet junction must declare which of its parameters
should be calibrated through a `calibrate` field listing the parameter names.
Parameters that are not listed are held constant at the value found in the
input file.

The names must match the parameter names a block exposes through its
`input_params` list (e.g. `R_poiseuille`, `C`, `L`, `stenosis_coefficient`
for `BloodVessel`).

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

### Resolution rules

* If a block has a `calibrate` field, only the listed parameters are
  optimized for that block.
* If a block has no `calibrate` field (or an empty list), every parameter of
  that block is held at its input value.
* The calibrator errors out if no parameter ends up selected anywhere in the
  model.

To calibrate every parameter of a block (the previous default), list every
name explicitly, e.g. `"calibrate": ["R_poiseuille", "C", "L", "stenosis_coefficient"]`
for a `BloodVessel` or `"calibrate": ["R_poiseuille", "L", "stenosis_coefficient"]`
for a `BloodVesselJunction`.

### Worked example

The repository ships a small fixture derived from the `0104_0001` Vascular
Model Repository case at
`tests/cases/vmr/input/0104_0001_calibrate_R_only.json`. Every vessel and
multi-outlet junction in that file carries `"calibrate": ["R_poiseuille"]`,
the non-R parameters are fixed at the calibrated reference, and every R
value is zeroed. The calibrator recovers the reference R values to machine
precision and leaves the other parameters untouched. The matching test is
`tests/test_calibrator.py::test_calibration_R_only`.

## Block requirements

A block is calibratable as long as it implements `update_gradient` and
exposes its parameter names through the standard `Block::input_params`
field. The calibrator reads this metadata at runtime via
`Block::input_params` and `Block::input_params_list`, so adding a new
calibratable block does not require any changes to `calibrate.cpp`.
