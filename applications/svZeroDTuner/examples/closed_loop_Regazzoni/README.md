# closed_loop_Regazzoni example

For full svZeroDTuner usage and configuration guidance, see the svZeroDTuner guide on the docs site: <https://simvascular.github.io/svZeroDSolver/tuner.html>.

This example tunes a four-chamber closed-loop heart model (based on Regazzoni et al.) with
systemic and pulmonary circulations. Scalar hemodynamic targets are used — systemic arterial
pressure and left-ventricular ejection fraction.

All inputs and targets are in SI units:

- Pressure: Pa
- Flow: m³/s
- Resistance: Pa·s/m³
- Compliance: m³/Pa

## Model

- Four `LinearElastanceChamber` blocks (LV, RV, LA, RA) with `piecewise_cosine` activation
- Systemic circuit: aortic valve → systemic arteries (`AR_SYS`) → systemic veins (`VEN_SYS`)
- Pulmonary circuit: pulmonary valve → pulmonary arteries (`AR_PUL`) → pulmonary veins (`VEN_PUL`)

Tunable parameters include cardiac elastances (`LV.Emax`, `RV.Emax`), vascular resistances
(`AR_SYS.R_poiseuille`), and systemic arterial compliance (`AR_SYS.C`).

## Tuning configurations

| File | Description |
|---|---|
| `tuning_nelder_mead.yaml` | Nelder-Mead; scalar systemic pressure + LV EF targets |
| `tuning_differential_evolution.yaml` | Differential evolution; same targets |
| `tuning_time_series_target.yaml` | Time-series aortic pressure target (`targets/target_pressure_ar_sys.csv`) |
| `tuning_complex.yaml` | Larger parameter set with log/max scaling |

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

2. (Optional) Sensitivity screening to identify influential parameters:

   ```bash
   svzerodtuner sensitivity-analysis sensitivity.yaml
   ```

   Writes `sensitivity_results/` with correlation scores and a ranked bar plot.

3. Optimization:

   ```bash
   svzerodtuner optimize tuning_nelder_mead.yaml
   ```

## Example outputs

### target_comparison.csv

Each row reports how closely the optimized model matches one target:

| Name | Target | Simulated | Target range | % Error |
|---|---|---|---|---|
| Systemic arterial max pressure | 13 065 Pa | ~13 065 Pa | ±10% | ~0 |
| Systemic arterial min pressure | 7 066 Pa | ~7 066 Pa | ±10% | ~0 |
| LV ejection fraction | 0.50 | ~0.50 | ±5% | ~0 |

### optimization_history/

| File | Contents |
|---|---|
| `objective_history.png` | Objective vs. iteration — confirm the curve flattens before termination |
| `parameter_evolution.png` | Per-parameter trajectories — flag parameters that converge against a bound |
| `history.csv` | Full numerical record for post-processing |

### sensitivity_results/

When sensitivity analysis is run first, `sensitivity_results/` contains:

| File | Contents |
|---|---|
| `sensitivity_scores.csv` | First-order correlation scores for each parameter × QoI pair |
| `sensitivity_bar_<qoi>.png` | Ranked bar chart showing which parameters drive each quantity most |

High-scoring parameters are the best candidates to include in the optimization; low-scoring
parameters can be fixed at baseline values to reduce the search space.

### optimized_simulation/

| File | Contents |
|---|---|
| `pressures.png` | Pressure waveforms at key junctions over the final cardiac cycle |
| `flows.png` | Flow waveforms through key vessels |
| `volumes.png` | Chamber volume waveforms (LV, RV, LA, RA) |
| `optimized_results.csv` | Full time-series for every output variable |
| `optimized_summary.csv` | Min / max / mean / std for every output variable |
