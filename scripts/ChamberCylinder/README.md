# ChamberCylinder: genet23 reproduction and verification

Scripts and data behind the `ChamberCylinder` block, the thick-walled
cylindrical left ventricle of Genet et al. (2023), "Reduced left ventricular
dynamics modeling based on a cylindrical assumption", Int. J. Numer. Meth.
Biomed. Engng., e3711 (genet23).

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

## Reference solutions

- **Digitized genet23 Figs 5-9** (`genet23_fig{5..9}_digitized.csv`,
  `ChamberCylinder_figures.py`):
  - Pressure, volume and the PV loop are matched in all 15 cases without
    fitting. The mean end-systolic volume error is 0.6 mL, 0.1 mL for Fig 5.
  - The torsion panels are not a valid reference (see below).
- **Independent solutions of the same equations**
  (`ChamberCylinder_verification.py`, no data from genet23):
  - symmetry;
  - static inflation against a semi-analytical solution;
  - free torsional, axial and breathing vibration against closed forms;
  - a full beat against an independent closed-loop re-implementation
    (agreement to 0.06 mL and 0.014° of twist).

  Its fast subset runs in `tests/test_chamber_cylinder_verification.py`.

## The genet23 inconsistency

**Old inertia term in the torsion figures.** The twist curves of genet23 Figs 5-9
were computed with an older CardiacLab inertia term:
- It stored the velocity integrated over the cylinder axis Z (Eq A7 as printed).
- This discards the Z-dependent twist and axial momentum every step. It acts
  like a damper of strength 2ΔM/Δt that slows the untwisting and depends on the
  time step.
- The code was corrected in March 2022, but the figures were not regenerated.
  Fig 4 (the convergence study) does come from the corrected code.

| Twist at 550 ms | Value |
| --- | --- |
| Fig 5 | 13.3° |
| Corrected CardiacLab | 5.3° |
| This block | 5.8° |

The authors confirmed this explanation.

**Remaining 0.5° difference: twist inertia.** Removing inertia gives 5.1° in
this block. The authors report that inertia changes their twist by less than
0.1°.

**Smaller issues:**
- **Valve capacitance:** Table 1's C_valve = 9e-9 m³/Pa is a typo for 9e-12.
  It is not modeled.
- **Missing inputs:** genet23 does not state the venous pressure, the initial
  arterial pressure or the input functions.

`ChamberCylinder_provenance.tex` documents the source of every parameter and
the full analysis.

## Files

| File | Content |
| --- | --- |
| `ChamberCylinder_genet_exact.py` | the genet23 baseline configuration |
| `ChamberCylinder_figures.py` | runs the cases of Figs 5-9 and overlays them on the digitized traces |
| `genet23_fig{5..9}_digitized.csv` | numeric traces digitized from the published figures |
| `ChamberCylinder_verification.py` | verification against independent solutions |
| `ChamberCylinder_provenance.tex` | parameter provenance, formulation audit, torsion analysis |
| `ChamberCylinder.yaml` | the continuous equations the block discretizes (documentation) |

Only digitized numbers are stored here; the paper's figures, text and the
authors' code are not.

## Running

Build the Python module (see the main README), then from the repository root:

```bash
PYTHONPATH=build/python python scripts/ChamberCylinder/ChamberCylinder_figures.py
PYTHONPATH=build/python python scripts/ChamberCylinder/ChamberCylinder_verification.py T1 T2 T3 T4 plot
```

Both write to `scripts/ChamberCylinder/output/` (not tracked). Build the
document from this folder with `pdflatex ChamberCylinder_provenance.tex`
(twice). The torsion diagnostics in the document used temporary kernel switches
that are not part of the block. They are archived in commit `bbd53847`, on the
branch `claude/chamber-cylinder-diagnostics-archive`. Figures that need them
appear as placeholders.
