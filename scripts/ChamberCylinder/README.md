# ChamberCylinder: genet23 setup

Setup of the `ChamberCylinder` block for the thick-walled cylindrical left
ventricle of Genet et al. (2023), "Reduced left ventricular dynamics modeling
based on a cylindrical assumption", Int. J. Numer. Meth. Biomed. Engng., e3711
(genet23).

## Setup

`ChamberCylinder_genet_exact.py` builds the genet23 baseline from existing
svZeroDSolver blocks:

- **Ventricle:** `ChamberCylinder` with the Table 1 geometry, the genet23
  passive material, tissue viscosity and the BCS active stress. The geometry is
  derived from the defining volumes (65 mL cavity, 130 mL wall, aspect 1.5).
- **Valves:** `PiecewiseValve` blocks (Eq 35).
- **Windkessel:** two-stage C-R-C-R chain of `BloodVessel` blocks (Eq 36).
- **Inputs:** the activation ν(t), atrial pressure P_at(t) and Frank-Starling
  curve n0(e_c) come from the authors' input functions. The venous pressure
  (500 Pa) and the initial arterial pressure (47 mmHg) are CardiacLab values
  supplied by the authors. No parameter is fitted.
- **Numerics:** time step 1 ms. The result is the fifth beat, as in genet23
  Sec 3.1. That beat is not periodic yet: the periodic end-systolic volume is
  0.6 mL larger.

The regression case `tests/cases/chamber_cylinder_genet_circulation.json` is
this setup on a coarse mesh (3 elements through the wall, 2 beats).

## Reference solution

The reference is genet23 Fig 5 and the current CardiacLab code of the authors.
With 12 elements through the wall:

| | EDV (mL) | ESV (mL) | Peak pressure (kPa) | Twist at 550 ms |
| --- | --- | --- | --- | --- |
| CardiacLab (current) | 135.7 | 74.1 | 12.79 | 5.3° |
| This block | 135.7 | 74.0 | 12.79 | 5.8° |

Pressure, volume and the PV loop also match the parameter studies of Figs 6-9
(15 cases, mean end-systolic volume error 0.6 mL) without fitting. The torsion
panels of Figs 5-9 are not a valid reference (see below).

## The genet23 inconsistency

**Old inertia term in the torsion figures.** The twist curves of genet23 Figs 5-9
were computed with an older CardiacLab inertia term:
- It stored the velocity integrated over the cylinder axis Z (Eq A7 as printed).
- This discards the Z-dependent twist and axial momentum every step. It acts
  like a damper of strength 2ΔM/Δt that slows the untwisting and depends on the
  time step.
- The code was corrected in March 2022, but the figures were not regenerated.
  Fig 4 (the convergence study) does come from the corrected code.

As a result, Fig 5 shows 13.3° of twist at 550 ms where the corrected code gives
5.3°. The authors confirmed this explanation.

**Remaining 0.5° difference: twist inertia.** Removing inertia gives 5.1° in
this block. The authors report that inertia changes their twist by less than
0.1°.

**Smaller issues:**
- **Valve capacitance:** Table 1's C_valve = 9e-9 m³/Pa is a typo for 9e-12.
  It is not modeled.
- **Missing inputs:** genet23 does not state the venous pressure, the initial
  arterial pressure or the input functions.

## Files

| File | Content |
| --- | --- |
| `ChamberCylinder_genet_exact.py` | the genet23 baseline configuration |
| `ChamberCylinder.yaml` | the continuous equations the block discretizes (documentation) |

## Running

Build the Python module (see the main README), then from the repository root:

```python
import sys
sys.path.insert(0, "scripts/ChamberCylinder")
import ChamberCylinder_genet_exact as G
import svzerod

result = svzerod.simulate(G.build())
```

`G.build()` returns the configuration as a dictionary; its arguments set the
mesh (`ne`), the number of beats (`ncycle`) and the inputs.
