
<div align="center">
<h1>svZeroDSolver</h1>

[![Test Status](https://github.com/simvascular/svZeroDSolver/actions/workflows/test.yml/badge.svg)](https://github.com/simvascular/svZeroDSolver/actions)
[![codecov](https://codecov.io/gh/SimVascular/svZeroDSolver/graph/badge.svg?token=FQKC9L5I0W)](https://codecov.io/gh/SimVascular/svZeroDSolver)
[![Latest Release](https://img.shields.io/github/v/release/simvascular/svZeroDSolver?label=latest)](https://github.com/simvascular/svZeroDSolver/releases/latest)
![Platform](https://img.shields.io/badge/platform-macOS%20|%20Ubuntu-blue)

</div>

svZeroDSolver is a fast simulation tool for modeling the hemodynamics of
vascular networks using zero-dimensional (0D) lumped parameter models.
You can find more information under the following links:

* [**Documentation**](https://simvascular.github.io/svZeroDSolver)
* [**Developer Guide**](https://simvascular.github.io/svZeroDSolver/developer_guide.html)
* [**svZeroDTuner Guide**](https://simvascular.github.io/svZeroDSolver/tuner.html)
* [**Bug Reports**](https://github.com/simvascular/svZeroDSolver/issues)
* [**Forum**](https://github.com/simvascular/svZeroDSolver/discussions)
* [**About SimVascular**](https://simvascular.github.io)

## Impedance BC (Olufsen periodic memory)

`IMPEDANCE` is a boundary condition that enforces

`P = Pd + z0 * Q + sum_{m=1..Nk-1} z[m] * Q_lag[m]`

with one-cycle memory stored in a ring buffer and coupling-safe trial/accept
semantics.

Required JSON fields in `bc_values`:

* `z`: time-domain impedance kernel array.

Optional fields:

* `Pd`: distal/reference pressure (default `0.0`).
* `convolution_mode`: `exact` (default) or `truncated`.
* `num_kernel_terms`: required when `convolution_mode` is `truncated`.

Simulation timing requirements:

* Non-coupled configs must set
  `simulation_parameters.number_of_time_pts_per_cardiac_cycle = z.size() + 1`.
* 3D-coupled configs must keep `simulation_parameters.number_of_time_pts = 2`
  for the external interface step. The canonical coupled simulation parameters
  are sufficient; if
  `simulation_parameters.number_of_time_pts_per_cardiac_cycle` is omitted, the
  impedance cycle discretization is inferred from
  `simulation_parameters.cardiac_period / dt`.
* In all cases, `simulation_parameters.cardiac_period / dt` must match
  `z.size()`.

For `truncated`, runtime cost is `O(num_kernel_terms)` per accepted 0D step.
