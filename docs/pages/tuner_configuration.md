@page tuner_configuration svZeroDTuner Configuration Reference

[TOC]

# Optimization Config Schema

```yaml
model:
  config_file: "model.json"

parameters:
  - name: "LV.Emax"
    bounds: [1e7, 1e9]
    scaling: log

targets:
  - name: Systemic arterial max pressure
    type: scalar
    expression: np.max(pressure:AV:AR_SYS)
    target_value: 13065
    relative_bounds: 5%
    weight: 1.0

objective:
  norm: L1

optimization:
  algorithm: "differential_evolution"
  terminate_at_zero: true
  maxiter: 200

output:
  directory: "optimization_results"
  save_history: true
  save_plots: true
  save_final_config: true
```

Required top-level sections:

- `model`
- `parameters`
- `targets`
- `objective`
- `optimization`

# Sensitivity Config Schema

```yaml
model:
  config_file: "model.json"

parameters:
  - name: "LV.Emax"
    bounds: [1e8, 5e8]

quantities_of_interest:
  - name: Systemic arterial max pressure
    expression: np.max(pressure:AV:AR_SYS)

sensitivity:
  n_samples: 256

output:
  directory: "sensitivity_results"
  save_plots: true
```

Required top-level sections:

- `model`
- `parameters`
- `quantities_of_interest`

# Field-by-Field Reference

## `model`

- `config_file` (`str`, required): path to svZeroD model JSON; relative paths are resolved from the YAML file location.

## `parameters[]`

- `name` (`str`, required): parameter key in `Block.Parameter` form.
- `bounds` (`list[2]`, required): lower/upper bounds with `min < max`.
- `scaling` (`str`, optional): one of `identity`, `log`, `max`.

## `targets[]` (optimization only)

- `name` (`str`, required): unique target label.
- `type` (`str`, required): `scalar` or `time_series`.
- `expression` (`str`, required): expression over solver outputs.
- `target_value` (`float`, scalar target option): point target.
- `target_range` (`list[2]`, scalar/time-series option): explicit target range `[min, max]`.
- `target_file` (`str`, time-series option): CSV with `time,value` columns.
- `relative_bounds` (`percent` or `[min,max]`, optional): range around target.
- `weight` (`float`, optional, default `1.0`): target weighting.

Legacy alias `uncertainty` is still recognized but should be replaced by `relative_bounds`.

## `objective`

- `norm` (`str`, required): `L1` or `L2`.

## `optimization`

- `algorithm` (`str`, required): `differential_evolution` or `Nelder-Mead`.
- `terminate_at_zero` (`bool`, optional, default `true`): stop early if objective reaches zero.
- Additional fields: forwarded directly to SciPy optimizer API.

## `output`

Optimization defaults:

- `directory`: `optimization_results`
- `save_history`: `true`
- `save_plots`: `true`
- `save_final_config`: `true`

Sensitivity default output directory:

- `directory`: `sensitivity_results`

## `quantities_of_interest[]` (sensitivity only)

- `name` (`str`, required): QoI label.
- `expression` (`str`, required): expression over outputs (scalar semantics).

## `sensitivity`

- `n_samples` (`int`, optional, default `512`): quasi-random sample count for screening.

# Validation Rules and Common Config Errors

Frequent validation failures include:

- Missing required sections (`model`, `parameters`, `targets`, `objective`, `optimization`)
- Invalid parameter bounds (`min >= max`)
- `log` scaling with non-positive bounds
- Scalar target missing both `target_value` and `target_range`
- Time-series target missing `target_file`
- Missing `objective.norm`
- Unknown target `type`
- Unknown optimizer `algorithm`

When tuning fails, start with [Troubleshooting](@ref tuner_troubleshooting).
