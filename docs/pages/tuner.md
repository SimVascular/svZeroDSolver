@page tuner svZeroDTuner Guide

[TOC]

# About

svZeroDTuner is a Python module and command-line tool for calibrating svZeroDSolver 0D model parameters against target hemodynamic quantities. It supports:

- Parameter optimization from YAML configuration files
- Scalar and time-series targets using expression-based output extraction
- Sensitivity analysis for screening influential parameters

The implementation is available in the `applications/svZeroDTuner` folder.

# When to Use / When Not to Use

Use svZeroDTuner when you need to:

- Fit uncertain model parameters to measured pressure, flow, or volume data
- Enforce physiologic target ranges rather than strict point matching
- Rank parameter influence before deciding what to calibrate

Do not use svZeroDTuner as a replacement for:

- 0D model authoring (use [svZeroDGUI](@ref GUI))
- Post-processing and network inspection (use [svZeroDVisualization](@ref visualization))
- Fundamental model-structure changes (update the model itself first)

# Quickstart

## CLI workflow

Run optimization:

```bash
svzerodtuner optimize applications/svZeroDTuner/examples/right_heart_pa/tuning_differential_evolution.yaml
```

Run sensitivity analysis:

```bash
svzerodtuner sensitivity-analysis applications/svZeroDTuner/examples/closed_loop_Regazzoni/sensitivity.yaml
```

Alias commands are also supported:

- `svzerodtuner run <config.yaml>` (alias for `optimize`)
- `svzerodtuner sensitivity <config.yaml>` (alias for `sensitivity-analysis`)

## Python API workflow

```python
from svzerodtuner.sv0d_tuner import SV0DTuner

# Optimization from YAML
result = SV0DTuner("applications/svZeroDTuner/examples/right_heart_pa/tuning_nelder_mead.yaml").optimize()
print(result["success"], result["best_value"])
```

# Workflow

A typical svZeroDTuner workflow is:

1. Run a baseline simulation and inspect available outputs.
2. Define tunable parameters and bounds.
3. Define target quantities (scalar and/or time-series).
4. Choose objective norm and optimization algorithm.
5. Run optimization and inspect history/termination diagnostics.
6. Validate the optimized model against targets and physiology.
7. Visualize model outputs and network behavior.

See [Worked Examples](@ref tuner_examples) for end-to-end templates.

# Related Tools

- [svZeroDVisualization Guide](@ref visualization) for plotting and network inspection.
- [svZeroDGUI Guide](@ref GUI) for graphical model construction.
- [svZeroDTuner Concepts](@ref tuner_concepts)
- [svZeroDTuner Configuration Reference](@ref tuner_configuration)
- [svZeroDTuner API Reference](@ref tuner_api)
- [svZeroDTuner Troubleshooting](@ref tuner_troubleshooting)
