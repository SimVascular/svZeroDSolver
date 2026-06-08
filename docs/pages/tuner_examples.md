@page tuner_examples svZeroDTuner Worked Examples

[TOC]

# Example 1: Minimal Scalar-Target Tuning

Recommended starting point:

- `applications/svZeroDTuner/examples/closed_loop_Regazzoni`

## Workflow

1. Baseline run

```bash
cd applications/svZeroDTuner/examples/closed_loop_Regazzoni
python -c 'from main import run_baseline; run_baseline("model.json")'
```

2. Configure scalar targets in one of:

- `tuning_differential_evolution.yaml`
- `tuning_nelder_mead.yaml`

3. Run optimization with the CLI

```bash
svzerodtuner optimize applications/svZeroDTuner/examples/closed_loop_Regazzoni/tuning_differential_evolution.yaml
```

## Expected artifacts

- `optimization_history/history.csv`
- `optimization_history/objective_history.png`
- `target_comparison.png` and `target_comparison.csv`
- optimized model JSON in output directory

# Example 2: Multi-Outlet Pulmonary Tree Tuning

Pulmonary tree example:

- `applications/svZeroDTuner/examples/right_heart_pa`

## Workflow

1. Run baseline and inspect generated baseline outputs:

```bash
cd applications/svZeroDTuner/examples/right_heart_pa
python -c 'from main import run_baseline; run_baseline("model.json")'
```

2. Tune pulmonary pressures and RPA/LPA flow split using:
- `tuning_differential_evolution.yaml`, or
- `tuning_nelder_mead.yaml`
3. Validate branch flow split and PA pressure range in output plots and CSV.

Run with CLI:

```bash
svzerodtuner optimize applications/svZeroDTuner/examples/right_heart_pa/tuning_nelder_mead.yaml
```

# Optional: Time-Series Target Matching

Time-series target example:

- `applications/svZeroDTuner/examples/closed_loop_Regazzoni/tuning_time_series_target.yaml`

This configuration demonstrates:

- `type: time_series`
- `target_file` with `time,value` columns
- range-based matching via `relative_bounds`

# Validation Checklist

For each example:

- Confirm optimization termination reason and success flag.
- Confirm objective trend decreases in `objective_history.png`.
- Compare target and simulated outputs in `target_comparison.csv`.
- Check whether best parameters are near bounds (printed warnings).
- Cross-check final waveforms with [svZeroDVisualization](@ref visualization) when needed.
