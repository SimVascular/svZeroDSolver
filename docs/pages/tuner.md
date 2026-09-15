@page tuner svZeroDTuner Guide

[TOC]

# About

svZeroDTuner is a Python module and command-line tool for calibrating svZeroDSolver 0D model parameters against target hemodynamic quantities. It supports:

- Parameter optimization from YAML configuration files
- Scalar and time-series targets using expression-based output extraction
- Sensitivity analysis for screening influential parameters

The implementation is available in the `applications/svZeroDTuner` folder.

# Additional Resources

- [svZeroDTuner Concepts](@ref tuner_concepts)
- [svZeroDTuner Configuration Reference](@ref tuner_configuration)
- [svZeroDTuner API Reference](@ref tuner_api)
- [svZeroDTuner Troubleshooting](@ref tuner_troubleshooting)

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

# When to Run Sensitivity Analysis

Sensitivity analysis in svZeroDTuner is a screening step, not a replacement for optimization.
It is most useful after you have a stable baseline model and a first draft of your candidate
parameters and quantities of interest, but before you launch a large optimization.

Run sensitivity analysis when:

- You have many plausible tunable parameters and need to decide which ones to include.
- You expect non-identifiability, meaning several parameters may produce similar changes in the targets.
- Optimization is expensive and you want to reduce the search space first.
- You are adding new targets or changing parameter bounds and want to re-rank parameter influence.

In practice, the recommended order is:

1. Run a baseline simulation and confirm the model solves cleanly.
2. Define the quantities of interest you care about matching.
3. Select a broad candidate parameter set with physically reasonable bounds.
4. Run sensitivity analysis to rank which parameters most affect those quantities.
5. Keep the high-impact parameters free, and fix or defer the low-impact ones.
6. Run optimization on the reduced parameter set.

Sensitivity analysis is especially valuable for larger closed-loop problems, where tuning too many
parameters at once can make optimization slow or poorly identifiable. It is usually optional for
small problems with only a few clearly relevant parameters.

Do not run sensitivity analysis as the very first step on an unvalidated model. If the baseline
simulation is unstable, has unit issues, or the quantities of interest are not yet defined, the
sensitivity ranking will not be reliable. You should also rerun it after major changes to the
parameter list, target definitions, or parameter bounds.

# Workflow

A typical svZeroDTuner workflow is:

1. Run a baseline simulation and inspect available outputs.
2. Define tunable parameters and bounds.
3. Define target quantities (scalar and/or time-series).
4. Optionally run sensitivity analysis to rank parameter influence and reduce the search space.
5. Choose objective norm and optimization algorithm.
6. Run optimization and inspect history/termination diagnostics.
7. Validate the optimized model against targets and physiology.
8. Visualize model outputs and network behavior.

See [Worked Examples](@ref tuner_examples) for end-to-end templates.

# Examples

Three worked examples are provided in `applications/svZeroDTuner/examples/`. Each has its own README with a model description, quickstart, and a table of expected output files.

| Example | Focus | README |
|---|---|---|
| `right_heart_pa` | Scalar pressure and RPA/LPA flow-split targets; right heart + reduced-order pulmonary tree | [README](../../applications/svZeroDTuner/examples/right_heart_pa/README.md) |
| `closed_loop_Regazzoni` | Scalar systemic pressure and LV ejection fraction targets; four-chamber closed-loop | [README](../../applications/svZeroDTuner/examples/closed_loop_Regazzoni/README.md) |
| `closed_loop_Zingaro` | Time-series chamber volume targets alongside scalar pressure and flow targets; 21-parameter closed-loop | [README](../../applications/svZeroDTuner/examples/closed_loop_Zingaro/README.md) |

See [svZeroDTuner Worked Examples](@ref tuner_examples) for step-by-step CLI walkthroughs.

# Related Tools

- [svZeroDVisualization Guide](@ref visualization) for plotting and network inspection.
- [svZeroDGUI Guide](@ref GUI) for graphical model construction.
