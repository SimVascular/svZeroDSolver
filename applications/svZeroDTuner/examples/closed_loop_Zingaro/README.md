# closed_loop_Zingaro example

For full svZeroDTuner usage and configuration guidance, see the svZeroDTuner guide on the docs site: <https://simvascular.github.io/svZeroDSolver/tuner.html>.

This example tunes a four-chamber closed-loop model against a richer target set that includes
**time-series chamber volume waveforms** derived from patient imaging data (`targets/P003_chamber_volumes/`),
alongside scalar pressure and flow targets. It also demonstrates log-scaled and max-scaled
parameter search spaces, which are recommended when parameter bounds span multiple orders of
magnitude or include variables that can be zero.

All inputs and targets are in SI units:

- Pressure: Pa
- Flow: m^3/s
- Volume: m^3
- Resistance: Pa*s/m^3
- Compliance: m^3/Pa

## Model

- Four `LinearElastanceChamber` blocks (LV, RV, LA, RA) with `piecewise_cosine` activation
- Systemic circuit: `R_UPSTREAM_SYS` -> `AR_SYS` -> `VEN_SYS`
- Pulmonary circuit: `R_UPSTREAM_PUL` -> `AR_PUL` -> `VEN_PUL`

Twenty-one parameters are tunable, covering cardiac elastances (`Emax`, `Epass`, `Vrest`),
vascular resistances and compliances, and an initial-condition pressure - all with `log` or
`max` scaling in `tuning.yaml`.

## Tuning configuration

`tuning.yaml` uses Nelder-Mead with `adaptive: true` and `maxiter: 10000`. The full parameter
list with recommended bounds and scaling is provided there; start with the smaller Regazzoni
example if you are new to svZeroDTuner.

## Quick start

Install the package once from the repository root:

```bash
pip install -e .
```

Then run from this directory:

1. Baseline inspection:

   ```bash
   python -c 'from main import run_baseline; run_baseline("model.json")'
   ```

   Writes `baseline_results/` with time-series CSV, summary CSV, and plots.

2. (Optional) Sensitivity screening:

   ```bash
   svzerodtuner sensitivity-analysis sensitivity.yaml
   ```

   Writes `sensitivity_results/` with correlation scores to help prioritize the 21 parameters.

3. Optimization:

   ```bash
   svzerodtuner optimize tuning.yaml
   ```

   This run is long (potentially hours) due to the number of parameters and time-series target
   evaluation cost. Use `workers: -1` with `differential_evolution` for parallel evaluation.

## Example outputs

### target_comparison.csv

For scalar targets, each row shows the target, simulated value, acceptable range, and percent
error. For time-series targets, a single summary row captures the mean absolute error across
all time points:

| Name | Type | Target | % Error |
|---|---|---|---|
| LV volume | time_series | patient waveform | <=5% at convergence |
| RV volume | time_series | patient waveform | <=5% at convergence |
| LA volume | time_series | patient waveform | <=5% at convergence |
| RA volume | time_series | patient waveform | <=5% at convergence |
| Aortic max pressure | scalar | 13 065 Pa | <=5% |
| Aortic min pressure | scalar | 7 066 Pa | <=5% |
| Aortic valve peak flow | scalar | 4.27 x 10^-4 m^3/s | <=20% |
| Pulmonary valve peak flow | scalar | 4.27 x 10^-4 m^3/s | <=20% |
| PA max pressure | scalar | 2 666 Pa | <=20% |
| PA min pressure | scalar | 1 533 Pa | <=20% |
| Systemic venous mean pressure | scalar | 800 Pa | <=20% |

### optimization_history/

| File | Contents |
|---|---|
| `objective_history.png` | Objective vs. iteration - for a 21-parameter run expect slow early progress |
| `parameter_evolution.png` | Per-parameter trajectories across all evaluations |
| `history.csv` | Full numerical record; useful for diagnosing stagnation or parameter bounds issues |

### sensitivity_results/

| File | Contents |
|---|---|
| `sensitivity_scores.csv` | First-order correlation scores for each of the 21 parameters x each QoI |
| `sensitivity_bar_<qoi>.png` | Ranked bar chart - use to identify which 5-8 parameters to tune first |

Running sensitivity analysis before full optimization is especially valuable with 21 parameters:
it can reduce a multi-hour run to a few minutes by fixing low-influence parameters.

### optimized_simulation/

| File | Contents |
|---|---|
| `pressures.png` | Pressure waveforms at key junctions over the final cardiac cycle |
| `flows.png` | Flow waveforms through key vessels (aortic valve, pulmonary valve) |
| `volumes.png` | Chamber volume waveforms overlaid with patient target curves |
| `optimized_results.csv` | Full time-series for every output variable |
| `optimized_summary.csv` | Min / max / mean / std for every output variable |
