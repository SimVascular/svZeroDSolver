# Valve convergence and robustness

Tools used to diagnose why valve blocks are expensive and fragile to solve,
and to compare the available valve formulations. Run everything from the
repository root with the compiled Python module on `PYTHONPATH`:

```bash
PYTHONPATH=$PWD/build/python python scripts/valve_convergence/<script>.py
```

| script | what it does |
| --- | --- |
| `diagnose.py` | Runs every test case that contains a valve and plots pressures across each valve, the valve state, chamber volume and the Newton iterations per time step. |
| `model.py` | Builds the physiologic benchmark: left atrial pressure source, mitral valve, elastance ventricle, aortic valve, three-element Windkessel. Parameterized by heart rate, contractility, pre- and afterload, valve type and unit system. |
| `benchmark.py` | Sweeps valve parameters, operating point, unit system, tolerance, line search and time-integrator damping, reporting both solver cost and physiologic accuracy. |
| `compare.py` | Figures comparing `ValveTanh`, `PiecewiseValve` and `ValveDiode` on the benchmark. |

The solver exposes the per-step convergence history used by these scripts:

```python
solver = svzerod.Solver(config)
solver.set_record_residual_history(True)   # optional, per-iteration residuals
solver.run()
solver.get_step_times()              # time at the start of each step
solver.get_nonlinear_iterations()    # Newton iterations per step
solver.get_step_residuals()          # largest residual left at the end
solver.get_step_residual_equations() # which block that residual belongs to
solver.get_step_residual_indices()   # and which equation index
solver.get_residual_history()        # residual after every Newton iteration
```
