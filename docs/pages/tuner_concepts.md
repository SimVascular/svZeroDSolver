@page tuner_concepts svZeroDTuner Concepts

[TOC]

# What Is Being Tuned

svZeroDTuner adjusts model parameters identified by name using the format `Block.Parameter`, for example:

- `LV.Emax`
- `AR_SYS.C`
- `RPA.R_poiseuille`

Names are resolved against svZeroD model JSON structures (chambers, vessels, valves, boundary conditions, and selected global sections).

# Targets and Metrics

Targets are defined in the YAML `targets` section and are computed from simulation outputs using expressions.

Two target types are supported:

- `scalar`: one value per target (for example `np.max(pressure:AV:AR_SYS)`)
- `time_series`: waveform target from a CSV file (`time`, `value`)

Expressions can combine outputs with `numpy` operations, and output names must match available solver output labels.

# Loss / Objective

The objective is built from weighted relative errors across all targets.

- `L1`: sum of absolute relative errors
- `L2`: Euclidean norm of the relative-error vector

Targets are internally treated as allowed ranges `[lo, hi]`:

- Point target: `target_value` implies `lo = hi = value`
- Relative range: `relative_bounds` expands around target value
- Explicit range: `target_range` provides `[min, max]`

Penalty is zero when simulated values are within range, and positive only when values are outside bounds.

# Constraints / Bounds / Scaling

Each parameter requires `bounds: [min, max]` and may define `scaling`:

- `identity`: no transform
- `log`: optimizer runs in log-space (bounds must be positive)
- `max`: scale by max bound magnitude

Bounds act as hard optimizer constraints and are also checked against initial parameter values.

# Convergence and Termination

svZeroDTuner currently supports:

- `differential_evolution`
- `Nelder-Mead`

Optimization options are passed through to SciPy using native option names.

By default, optimization can terminate early when objective reaches zero-range penalty (`terminate_at_zero: true`).

# Failure Modes

Common failure patterns:

- Simulation failure for trial parameters (objective penalized with a large fallback value)
- Multiple parameter sets producing similar outputs (non-identifiability)
- Persistent convergence at parameter bounds (model or bounds may be restrictive)
- Conflicting targets that cannot be satisfied simultaneously

See [Troubleshooting](@ref tuner_troubleshooting) for mitigation strategies.
