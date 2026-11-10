# ChamberCylinder: verification against 3D finite elements

The 0D `ChamberCylinder` block is compared with svMultiPhysics solving the same
cylinder in 3D. Each test case is defined by a 0D input file; `cylinder_3d.py`
turns it into the matching 3D case, so the two always share geometry, material
and load.

| Step | Case | Status |
| --- | --- | --- |
| 1 | `passive_isotropic/`: quasi-static inflation, neo-Hookean wall | implemented |
| 2 | `passive_anisotropic/`: the same wall with one circumferential fiber family (Holzapfel-Ogden) | implemented |
| 3 | PV loop: 0D and 3D cylinders coupled to the same circulation | planned |

## Problem definition

- **Fixed length.** Both end faces keep u_z = 0 and slide freely in their
  planes; the 0D block uses `fixed_length = 1` (ε = 0, unloaded lids). The 0D
  default (lids loaded by the cavity pressure) has no 3D counterpart in
  svMultiPhysics: the end face would have to stay planar while carrying a force
  that depends on the unknown inner radius.
- **Eighth model.** The 0D solution is uniform along the axis and in θ, so the
  3D mesh covers a 90° sector with symmetry planes (u_y = 0 on θ = 0, u_x = 0
  on θ = 90°) and half the length (u_z = 0 on the mid-length plane). There are
  no rigid-body modes and no stabilizing springs.
- **Fiber architecture.** The symmetry planes and the free-sliding ends need
  σ_θZ = 0, which holds for circumferential or axial fibers. Helical fibers
  make the cylinder twist and load every cross-section with a through-wall
  shear σ_θZ(R) of zero resultant; the uniform 0D solution then needs ends that
  rotate rigidly while expanding radially (genet23's lids). svMultiPhysics has
  no such boundary condition (`Undeforming_neu_face` forces one displacement
  for the whole face), so a full 3D cylinder with helical fibers develops end
  layers and does not converge to the 0D model. `cylinder_fibers.py` writes
  the 0D fiber and sheet (e_R) directions per node; svMultiPhysics
  interpolates them to the quadrature points, so they follow the curved basis
  vectors to the order of the elements. Per-element directions
  (`fiber_location="cells"`) are constant in each element, an O(h) error that
  limits the 3D convergence.
- **Materials.** The 3D side uses svMultiPhysics' material laws unchanged;
  where the 0D law differs, the 0D side adopts the 3D law. With C3 = 0, the 0D
  passive law is the svMultiPhysics neo-Hookean (C1 = C10 = μ/2), Mooney-Rivlin
  (C1, C2) or HGO fiber term (a4 = 2 C5 C6, b4 = C6, κ = 0). The `ho_*`
  parameters add svMultiPhysics' `HolzapfelOgden` law (isotropic, fiber, sheet
  and fiber-sheet terms, sheet = e_R, smoothed Heaviside slope k).
  `check_materials.py` compares the two stress implementations for random
  deformations (agreement to round-off). Viscosity γ = 0.
- **Volumetric law.** The 0D input sets `bulk_modulus` (κ of svMultiPhysics'
  ST91 energy U = κ/4 (J² − 1 − 2 ln J)) and both models use it, with
  κ = 100× the isotropic shear modulus as in the svMultiPhysics struct tests.
- **Quasi-static.** The 0D runs without inertia (`density` = 0) and the 3D with
  a negligible density. The pressure ramp is slow because the 0D active branch
  acts as a fiber dashpot μ even without activation; a 1200 s ramp makes its
  stress ~1e-6 of the elastic stress.
- **Load and units.** The cavity pressure spans 0–120 mmHg (end-diastolic to
  peak-systolic levels) in 24 steps of 50 s (5 mmHg). Solver inputs are SI
  (Pa); reports and figures use mmHg (1 mmHg = 133.322 Pa).

### Anisotropic case

The wall of step 1 (μ = 100 kPa) with one circumferential fiber family,
`HolzapfelOgden` with a = 100 kPa, b = 0 (neo-Hookean), a_f = 100 kPa,
b_f = 1, k = 0 and the sheet terms off. Two choices follow from tests of the
svMultiPhysics default parameters (struct/LV_HolzapfelOgden_passive):

- k = 0. With k = 100 the smoothed Heaviside makes the fiber and sheet
  energies non-convex just below zero strain (the sheet term's stiffness drops
  to −400 Pa at −1.6% strain, 7× the isotropic a = 59 Pa). The 3D problem then
  loses ellipticity in the first mmHg of inflation and Newton diverges.
- Mild anisotropy. The default a = 59 Pa is 300× softer than the fibers, so a
  penalty of 100 a leaves the wall highly compressible (wall volume ×1.9 at
  120 mmHg), while a penalty of 100 a_f makes the linear systems too stiff for
  the fsils solvers.

The 0D Holzapfel-Ogden implementation is also checked against the closed
form of the incompressible tube (`fixed_length_exact_ho`,
`tests/test_chamber_cylinder_verification.py`).

## Compared quantities

`run_study.py` compares, at every 0D time point of the ramp, the 0D kinematic
quantities with the 3D solution, averaged over the nodes at each radial node
level R_k (the levels are nodes of the 0D reference mesh):

| Key | Quantity | Error measure |
| --- | --- | --- |
| `V` | cavity volume | max \|V − V_ref\| / V_ref |
| `rho` | radial displacement ρ(R) | max \|ρ − ρ_ref\| / max \|ρ_ref\| |
| `h` | wall-thickness change ρ(Re) − ρ(Ri) | relative to its maximum |
| `beta` | twist per unit length β (helical fibers) | relative to max \|β_ref\| |
| `phi` | transmural rotation φ(R) − φ(Ri) | relative to max \|β_ref\| L |
| `eta` | axial displacement η(R) − η(Ri) | relative to max \|ρ_ref\| |

The reference is the closed form for the incompressible isotropic tube and
otherwise the Richardson extrapolation of the 0D block at `num_elements` and
twice that (the 0D model converges at second order). The study also reports
the 0D convergence against the same reference and the non-homogeneity of the
3D solution (largest deviation of ρ from its level mean).

### Struct output fields

svMultiPhysics writes Displacement, Velocity, Acceleration and
Fiber_direction. Its other struct outputs (Stress, Cauchy_stress,
VonMises_stress, Strain, Def_grad, Jacobian, Fiber_alignment, Volume) abort
for quadratic elements with mid-edge nodes on curved surfaces (`get_nnx`:
"Error in computing shape functions"), so `cylinder_3d.py post` computes them
from the displacement with the curved element geometry, averaged to the nodes:
deformation gradient, Jacobian, Green-Lagrange strain, 2nd Piola-Kirchhoff and
Cauchy stress (the svMultiPhysics material model of solver.xml), von Mises
stress, fiber stretch and deformed fiber direction, and the cylindrical
components of strain and Cauchy stress. They are written per step to
`CASE/runs/<run>/post/fields_<step>.vtu` (ParaView), and their through-wall
profiles (J, E_RR, E_TT, E_ZZ, σ_rr, σ_tt, σ_zz, σ_vm, λ_f) are in
`kinematics.json`. With κ = 100 μ the nodal hydrostatic stress oscillates
between node levels (about ±2 kPa at nr = 4): it is κ (J − 1) of a nodal J
from discontinuous gradients. Strains, fiber stretch and von Mises stress are
smooth.

## 3D formulations (`--formulation`)

| Name | Elements | Volumetric law | Notes |
| --- | --- | --- | --- |
| `struct-p2` (default) | P2 displacement | ST91 penalty κ = 0D `bulk_modulus` | |
| `struct-p1` | P1 displacement | ST91 penalty | locks as κ/μ grows |
| `ustruct-th` | P2-P1 Taylor-Hood | incompressible (ν = 0.5) | no svMultiPhysics test case uses it |
| `ustruct-p1` | P1-P1, VMS-stabilized | incompressible (ν = 0.5) | stabilization error ∝ h |

Solver settings follow the svMultiPhysics struct tests: Newton with 3 to 20
iterations and tolerance 1e-10 (LV_HolzapfelOgden_passive), BiCGStab with the
fsils preconditioner, tolerance 1e-12 and at most 600 iterations
(block_compression). The linear absolute tolerance is 1e-20 because SI
residuals are small. Every run streams the solver output and stops at once on
NaN, a diverging or unconverged Newton iteration, or a time limit
(`--max-time`); `--substeps` refines the 3D load steps.

## Running

Build the Python module (main README), then from the repository root:

```bash
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/run_study.py scripts/ChamberCylinder/verification_3d/passive_anisotropic --nr 1 2 4 --run --np 4
```

`--nr` sets the elements through the wall (the mesh-size parameter;
circumferential and axial divisions follow from h = (Re − Ri)/nr). Cases go to
`CASE/runs/`, results and figures to `CASE/output/`, which also caches the 0D
runs. Without `--run` the script only generates the 3D cases and
post-processes existing results.

The 3D solver runs in the `simvascular/solver:latest` Docker image (or a
local build with `--exe`). The image is built for x86-64; on Apple Silicon,
Docker runs it under emulation. Single tools:

```bash
python scripts/ChamberCylinder/verification_3d/cylinder_3d.py generate MODEL_0D.json CASE_DIR --nr 4 --half-length
python scripts/ChamberCylinder/verification_3d/cylinder_3d.py run CASE_DIR --np 4
python scripts/ChamberCylinder/verification_3d/cylinder_3d.py post CASE_DIR
python scripts/ChamberCylinder/verification_3d/check_materials.py MODEL_0D.json CASE_DIR
```

Post-processing needs VTK (`pip install vtk`).
