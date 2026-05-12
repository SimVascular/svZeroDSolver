@page tuner_troubleshooting svZeroDTuner Troubleshooting

[TOC]

# Numerical Instability

Symptoms:

- Frequent simulation failures during optimization
- Objective stuck at a very large fallback value

Actions:

- Narrow parameter bounds to physically plausible ranges.
- Switch to conservative optimizer settings (for example Nelder-Mead with tighter search region).
- Reduce number of simultaneously tuned parameters.

# Non-Identifiability / Over-Parameterization

Symptoms:

- Many parameter sets yield similar objective values
- Parameters drift to bounds with little objective change

Actions:

- Add more independent targets.
- Fix insensitive parameters using sensitivity analysis.
- Tune in stages (global resistances/compliances first, then chamber or valve parameters).

# Bad Initial Guess or Bounds

Symptoms:

- Immediate validation error for out-of-bounds initial values
- Early optimizer stagnation

Actions:

- Check model JSON initial values against tuning bounds.
- Expand or shift bounds where justified.
- Use scaling (`log`/`max`) for wide-dynamic-range parameters.

# Inconsistent or Conflicting Targets

Symptoms:

- Persistent nonzero objective despite long runs
- One target improves while another worsens

Actions:

- Re-check units and physiological consistency.
- Use target ranges (`relative_bounds` or `target_range`) instead of exact points.
- Rebalance target `weight` values.

# Coupling / Units Mismatch

Symptoms:

- Unphysical magnitudes (for example pressure or flow by orders of magnitude)

Actions:

- Verify model and target units are consistent (SI vs cgs).
- Confirm expression output names and expected dimensions.
- Validate baseline outputs before tuning.

# Expression Errors

Symptoms:

- Expression evaluation exceptions
- Missing output name messages

Actions:

- Use exact solver output names from baseline result summaries.
- Keep expressions limited to available outputs and valid `numpy` syntax.
- For time-series targets, verify CSV has `time` and `value` columns.

# Interpreting Bound Warnings

svZeroDTuner warns when best parameters are close to bounds.

Interpretation:

- The optimum may be constrained by bounds.
- The model may require wider bounds or a different parameterization.

Actions:

- Expand bounds only when physically justified.
- Add constraints via additional targets or tighter physiologic ranges.

# Debug Checklist

1. Validate config sections and required fields.
2. Run baseline first and confirm output names.
3. Verify target file format (`time,value`) and units.
4. Start with fewer parameters and scalar targets.
5. Inspect `history.csv` and `target_comparison.csv` after each run.
6. Run sensitivity analysis to identify low-impact parameters.
