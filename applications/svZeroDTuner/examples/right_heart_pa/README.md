# right_heart_pa example

For full svZeroDTuner usage and configuration guidance, see the svZeroDTuner guide on the docs site: <https://simvascular.github.io/svZeroDSolver/tuner.html>.

This example couples a prescribed venous inflow waveform to a right atrium and right ventricle model using `LinearElastanceChamber` blocks with `piecewise_cosine` activation and `PiecewiseValve` blocks. The right ventricle ejects into a reduced-order pulmonary artery model with:

- MPA: constant resistance
- RPA/LPA: resistance + stenosis coefficient + inductance
- RPA/LPA outlets: RCR boundary conditions
- Upstream venous return: `VEN_SYS` vessel (CRL) driven by `VEN_IN_FLOW`

All inputs and targets in this example are in cgs units:

- Pressure: dyn/cm^2
- Flow: cm^3/s
- Resistance: dyn\*s/cm^5
- Compliance: cm^5/dyn
- Inductance: dyn\*s^2/cm^5

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

Install the package once from the repository root so the `svzerodtuner` command is available:

```bash
pip install -e .
```

Then run the example from this directory:

1. Baseline inspection:

   ```bash
   python -c 'from main import run_baseline; run_baseline("model.json")'
   ```

   This creates `baseline_results/` with the full time series, summary CSV, and plots for target selection.

2. Tuning:
   - Update targets in the tuning YAML if needed.
   - Run one of the tuning configurations with the CLI:

   ```bash
   svzerodtuner optimize tuning_nelder_mead.yaml
   ```

   or

   ```bash
   svzerodtuner optimize tuning_differential_evolution.yaml
   ```

   `svzerodtuner run <config.yaml>` is also supported as an alias for `optimize`.

## Notes

Target names follow the solver output naming scheme:

- `pressure:<block1>:<block2>`
- `flow:<block1>:<block2>`

Cardiac output is defined as mean flow through `PV -> MPA`.
Flow split is enforced by targeting the mean flows through `RPA -> RCR_RPA` and `LPA -> RCR_LPA`.

## Example outputs

After running `svzerodtuner optimize tuning_nelder_mead.yaml`, the output directory `optimization_results_right_heart_pa_nm/` contains the following.

### target_comparison.csv

Each row reports the target, simulated value, acceptable range, and percent error for one quantity:

| Name | Target | Simulated | Target range | % Error |
|---|---|---|---|---|
| Pulmonary arterial mean flow | 83.3 cm³/s | 84.3 cm³/s | [75.0, 91.7] | 1.2 |
| Pulmonary arterial max pressure | 33 330 dyn/cm² | 36 660 dyn/cm² | [29 997, 36 663] | 10.0 |
| Pulmonary arterial min pressure | 10 670 dyn/cm² | 11 728 dyn/cm² | [9 603, 11 737] | 9.9 |
| Pulmonary arterial mean pressure | 20 000 dyn/cm² | 18 000 dyn/cm² | [18 000, 22 000] | -10.0 |
| RPA mean flow | 45.8 cm³/s | 42.3 cm³/s | [41.2, 50.4] | -7.8 |
| LPA mean flow | 37.5 cm³/s | 38.3 cm³/s | [33.8, 41.3] | 2.0 |

All targets land within their acceptable ranges.

### final_parameters.json

```json
{
  "RA.Emax": 91.5,   "RA.Epass": 58.2,
  "RV.Emax": 358.0,  "RV.Epass": 60.6,
  "MPA.R_poiseuille": 31.0,
  "RPA.R_poiseuille": 61.4,  "RPA.L": 0.607,  "RPA.stenosis_coefficient": 3.1e-4,
  "LPA.R_poiseuille": 70.3,  "LPA.L": 0.450,  "LPA.stenosis_coefficient": 2.6e-4,
  "RCR_RPA.Rp": 161.1,  "RCR_RPA.Rd": 146.2,  "RCR_RPA.C": 7.3e-3,
  "RCR_LPA.Rp": 177.8,  "RCR_LPA.Rd": 158.8,  "RCR_LPA.C": 6.6e-3
}
```

### optimization_history/

| File | Contents |
|---|---|
| `history.csv` | Objective value and all parameter values at every function evaluation |
| `objective_history.png` | Objective vs. iteration — use to confirm convergence |
| `parameter_evolution.png` | Each parameter value vs. iteration — use to spot parameters hitting bounds |
| `timing_info.json` | Total wall time and per-evaluation average (789 evals, ~4 min 39 s for this run) |

### optimized_simulation/

| File | Contents |
|---|---|
| `pressures.png` | Pressure waveforms at key junctions over the final cardiac cycle |
| `flows.png` | Flow waveforms through key vessels |
| `volumes.png` | Chamber volume waveforms (RA, RV) |
| `optimized_results.csv` | Full time-series for every output variable |
| `optimized_summary.csv` | Min / max / mean / std for every output variable |
