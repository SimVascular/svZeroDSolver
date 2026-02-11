# right_heart_pa example

This example couples a prescribed venous inflow waveform to a right atrium and right ventricle model using `PiecewiseCosineChamber` and `PiecewiseValve` blocks. The right ventricle ejects into a reduced-order pulmonary artery model with:

- MPA: constant resistance
- RPA/LPA: resistance + stenosis coefficient + inductance
- RPA/LPA outlets: RCR boundary conditions
- Upstream venous return: `VEN_SYS` vessel (CRL) driven by `VEN_IN_FLOW`

All inputs and targets in this example are in cgs units:
- Pressure: dyn/cm^2
- Flow: cm^3/s
- Resistance: dyn*s/cm^5
- Compliance: cm^5/dyn
- Inductance: dyn*s^2/cm^5

## Current Model Notes

- The inflow BC `VEN_IN_FLOW` is populated from the Regazzoni closed-loop baseline `flow:VEN_SYS:J1` waveform (converted to cgs).
- IVC/SVC have been replaced by a single venous inflow that feeds `VEN_SYS` and `J_VEN -> RA`.
- Cardiac period is set to **0.689 s** to match the Regazzoni waveform.
- Pulmonary compliance has been increased relative to the initial MPA/RPA/LPA settings.
- Solver settings are more conservative (higher nonlinear iterations, smaller time step).

Two tuning configs are provided:
- `tuning_differential_evolution.yaml`
- `tuning_nelder_mead.yaml`

In practice, Nelder-Mead has produced more stable convergence for this model than differential evolution.

## Quick start

1. Baseline:
   - Edit `main.py` and uncomment `run_baseline("model.json")`.
   - Run: `python main.py`

2. Tuning:
   - Update targets in the tuning yaml if needed.
   - Uncomment one of the `run_optimization(...)` lines in `main.py`.
   - Run: `python main.py`

## Notes

Target names follow the solver output naming scheme:
- `pressure:<block1>:<block2>`
- `flow:<block1>:<block2>`

Cardiac output is defined as mean flow through `PV -> MPA`.
Flow split is enforced by targeting the mean flows through `RPA -> RCR_RPA` and `LPA -> RCR_LPA`.
