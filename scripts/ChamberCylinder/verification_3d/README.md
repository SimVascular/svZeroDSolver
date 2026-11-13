# ChamberCylinder: verification against 3D finite elements

The 0D `ChamberCylinder` block is compared with svMultiPhysics solving the same
cylinder in 3D. Each test case is defined by a 0D input file; `cylinder_3d.py`
turns it into the matching 3D case, so the two always share geometry, material
and load.

| Step | Case | Status |
| --- | --- | --- |
| 1 | `passive_isotropic/`: quasi-static inflation, neo-Hookean wall | implemented |
| 2 | `passive_anisotropic/`: the same wall with one circumferential fiber family (Holzapfel-Ogden) | implemented |
| 2b | `passive_anisotropic_full/`: all Holzapfel-Ogden terms (exponential isotropic, fiber, sheet, fiber-sheet) with the tension-compression switch | implemented |
| 2b' | `passive_isotropic_closed/`, `passive_anisotropic_closed/`, `passive_anisotropic_full_closed/`: the walls of steps 1–2b with free length (0D ends `lids`, loaded by the cavity pressure), 3D tube closed by an end plate, 6 L0 long | implemented (nr 1, 2, 4; full HO 2, 4, 8) |
| 2c | `passive_helical_L2/`: the wall of step 2 with helical fibers (+60° to −60°), full cylinder with free ends, twice the test-case length; `passive_helical_a30_L2/` (uniform 30°) and `passive_helical_a60to0_L2/` (+60° to 0°) | implemented (nr 1, 2, 4) |
| 3a | `active_prescribed/`: PV loop with a prescribed active tension, 3D cylinder coupled to the circulation of `tests/cases/chamber_cylinder.json` | implemented (nr 1, 2, 4, 8) |
| 3a' | `active_helical_L2/`: the same with helical fibers (+60° to −60°), full 3D cylinder with free ends, twice the length | implemented (nr 1, 2) |
| 3b | PV loop with an active-stress ODE model (e.g. Land-Niederer in both) | planned |

## Problem definition

- **0D formulation.** The block's `formulation` object
  (`src/model/ChamberCylinder.h`) swaps each component between genet23 (the
  default) and its svMultiPhysics equivalent: `passive` (`genet23` C1–C6 or
  `holzapfel_ogden`), `volumetric` (`incompressible` or `st91`), `active`
  (`bcs` or `uniform_unsteady`) and `ends` (`lids` or `fixed_length`). Every
  case input states all four. The passive cases keep the BCS law at rest
  (`sigma_max` = 0).
- **Length.** Except in the free-length cases, both end faces keep u_z = 0
  and slide freely in their planes; the 0D block uses the ends
  `fixed_length` (ε = 0, unloaded lids). The free-length cases (`*_closed`)
  use the genet23 ends `lids`, loaded by the cavity pressure, and close the
  3D tube with a meshed end plate (see "Free length" below).
- **Eighth model.** The 0D solution is uniform along the axis and in θ, so the
  3D mesh covers a 90° sector with symmetry planes (u_y = 0 on θ = 0, u_x = 0
  on θ = 90°) and half the length (u_z = 0 on the mid-length plane). There are
  no rigid-body modes and no stabilizing springs.
- **Fiber architecture.** The symmetry planes and the free-sliding ends need
  σ_θZ = 0, which holds for circumferential or axial fibers. Helical fibers
  make the cylinder twist and load every cross-section with a through-wall
  shear σ_θZ(R) of zero resultant; the 3D model then covers the full cylinder
  with free ends and differs from the 0D near the ends (see "Helical fibers"
  below). `cylinder_fibers.py` writes
  the 0D fiber and sheet (e_R) directions per node; svMultiPhysics
  interpolates them to the quadrature points, so they follow the curved basis
  vectors to the order of the elements. Per-element directions
  (`fiber_location="cells"`) are constant in each element, an O(h) error that
  limits the 3D convergence.
- **Materials.** The 3D side uses svMultiPhysics' material laws unchanged;
  where the 0D law differs, the 0D side adopts the 3D law. With C3 = 0, the 0D
  passive law is the svMultiPhysics neo-Hookean (C1 = C10 = μ/2), Mooney-Rivlin
  (C1, C2) or HGO fiber term (a4 = 2 C5 C6, b4 = C6, κ = 0). The passive law
  `holzapfel_ogden` is svMultiPhysics' `HolzapfelOgden` law with the `ho_*`
  parameters (isotropic, fiber, sheet and fiber-sheet terms, sheet = e_R,
  smoothed Heaviside slope k).
  `check_materials.py` compares the two stress implementations for random
  deformations (agreement to round-off). Viscosity γ = 0.
- **Volumetric law.** The 0D input selects the volumetric law `st91` with
  `bulk_modulus` (κ of svMultiPhysics' ST91 energy
  U = κ/4 (J² − 1 − 2 ln J)) and both models use it, with
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

`passive_anisotropic_full/` uses every term of the law: a = 50 kPa, b = 5,
a_f = 25 kPa, b_f = 20, a_s = 25 kPa, b_s = 5, a_fs = 10 kPa, b_fs = 5, the
svMultiPhysics default k = 100 and κ = 100 a = 5 MPa. Inflation stretches the
circumferential fibers and compresses the radial sheets, so the switch turns
the sheet term off; the fiber exponential halves the compliance over the ramp
(+26% lumen volume at 120 mmHg). With these ratios the non-convex dip of the
switched terms (−4 kPa at −1.6% strain) stays far below the isotropic
stiffness. The fiber-sheet term carries no stress in the exact solution
(I8 = 0) but enters the tangent. Run it with `--nr 2 4` or finer: with one
element through the wall (nr = 1) the discrete solution leaves the uniform
branch above 70 mmHg (non-homogeneity 4% at 70 mmHg, 10% at 75 mmHg, 990%
at 120 mmHg), an
instability of the coarse mesh triggered by the switch (with k = 0 the same
mesh stays uniform), while nr = 2 and 4 stay uniform to within 1%.

The 0D Holzapfel-Ogden implementation is also checked against the closed
form of the incompressible tube (`fixed_length_exact_ho`,
`tests/test_chamber_cylinder_verification.py`).

### Helical fibers

`passive_helical_L2/` is the wall of `passive_anisotropic/` with the helix of
`tests/cases/chamber_cylinder.json`, +60° at the endocardium to −60° at the
epicardium, in a tube twice as long as the test case (2 L0 = 114 mm, initial
volume doubled). Mirroring θ or Z flips the helix angle, so the eighth
model's symmetry planes do not apply: the 3D model is the full cylinder over
the full length. It is as free as the 0D: both ends keep u_z = 0 (fixed
length) and are otherwise free to twist and expand. In-plane dashpots on both
end faces (Robin BCs with zero stiffness, damping 1e-3 μ/Re·dt) hold the three
rigid-body modes (x, y translation, rotation about z): within a load step they
act as springs of 1e-3 μ/Re, at rest they carry no force (during the ramp
< 2e-5 of the cavity pressure). No 0D data enters the 3D model. Alternatives
that did worse: a weak epicardial spring (1e-6 μ/Re) gives the same result but
175–950 linear iterations per solve; no support at all with GMRES stalls
(Newton stops at 4e-9, the discrete system is not exactly consistent); two
pinned radial lines at the base give the same averages but disturb the wall
near the pins.

The uniform 0D state does not solve this 3D problem exactly. Every 0D
cross-section carries the through-wall shear σ_θZ(R) (positive near the
endocardium, negative near the epicardium, zero net torque, up to 18% of
σ_θθ at 120 mmHg); the free 3D ends carry none. In end layers the 3D wall
twists faster and expands up to 1.6% more. `run_study.py` therefore compares
helical cases at the mid-length section only (the middle element layer).
`helical_sketch.py` draws both models and `helical_figure.py` the profiles
along the axis, with the 0D section shear and twist against the helix angle
(`output/helical_bcs_sketch.pdf`, `output/helical_discrepancy.pdf`).

Why twice the length: in the test-case tube (L0 = 57 mm, about four wall
thicknesses) the end layers reach the middle. At nr = 4 its mid-length
endocardial ρ differs from that of the 2 L0 tube by 3.4e-3 and its twist rate
by 14%, and neither converges to the 0D. A tube four times as long gives the
mid-length state of the 2 L0 tube to 1e-4 in ρ at nr = 2 (2e-5 for the uniform
30° helix), so the middle of the 2 L0 tube is free of the end effect. Ten times
or a tenth of the dashpot damping changes the mid-length ρ by 1.5e-5 and
1.5e-6, and dashpots on the epicardium instead of the ends give the same state
to 1e-5. (These checks are in `passive_helical/`, `passive_helical_L4/`,
`passive_helical_a30/`, `passive_helical_a30_L4/` and `active_helical/`.)

The full cylinder is expensive: over 2 L0 it has about 15× the elements of
the eighth model at the same nr, and it is held only weakly against rigid motion. The
fsils BiCGStab solver stalls from nr = 2 on, with any support; Trilinos GMRES
with ML (`--ls trilinos-ml`) converges with 500–1000 iterations per solve.
nr = 4 (38016 elements) takes 5 h on the cluster with 2 MPI ranks.

### Free length (closed ends)

With the ends `lids` (the genet23 default) the 0D block closes the tube with
lids that carry the cavity pressure and leaves the axial strain ε free
(genet23 Eqs. 39–44): the wall carries the axial force P π r_i² and
V = π r_i² (1 + ε) L. svMultiPhysics' caps (`Chamber_cap_surface`) close the
cavity for the flux but carry no pressure, so `generate(..., closed_end=True)`
meshes an end plate (`mesh_closed_cylinder`): the eighth model from the
mid-length plane (Z = 0, u_z = 0) to Z = L/2, continued by a plate of the wall
thickness H = Re − Ri that covers the cavity (domain 2 in `DOMAIN_ID`, an
isotropic neo-Hookean with the small-strain shear and bulk moduli of the
wall). The follower cavity pressure acts on the endocardium (`endo`) and the
inner face of the plate (`cap`); nothing holds the length. By equilibrium the
wall carries the axial force of the 0D lids, and no 0D data enters the 3D
model. `run_study.py` picks this setup when the 0D formulation has the ends
`lids`.

The plate holds the end ring in (endocardial ρ 80% below the 0D value at the
plate). This end layer decays with an oscillation: in the isotropic wall at
nr = 4 the ρ deviation is 7% at 2 H from the plate, 0.5% at 3 H, and at the
discretization level (1e-3) from about 6 H on; the axial strain, relative to
its small value (0.5% at 120 mmHg), is still 9% off at 3–4 H and within 1%
from 6 H on. The test-case tube (L0 = 4.1 H) would put the mid-length plane
inside this layer, so both models use a tube six times as long (6 L0 = 343
mm, 12.3 H from the plate to the mid-length plane, initial volume 6 V0). The
mid-length plane is compared: the first element layer, ε from a fit of u_z,
V = π (Ri + ρ(Ri))² (1 + ε) L as in the 0D. `closed_end_figure.py` draws the
end layer against the tube length and the mesh
(`passive_isotropic_closed/output/closed_end.pdf`).

### Active contraction in a circulation

`active_prescribed/` is the 0D test case `tests/cases/chamber_cylinder.json`
(pressure source, tanh valves, RC vessels, 1 s period) with the cylinder
adapted to what svMultiPhysics can represent: fixed length, circumferential
fibers, the Holzapfel-Ogden wall of `passive_anisotropic/`, no viscosity, no
wall inertia, and the BCS active law replaced by a prescribed tension
T(t) = 80 kPa sin²(π (t − 0.1 s)/0.4 s) on 0.1–0.5 s. The 0D active law
`uniform_unsteady` with its `active_tension` object reproduces
svMultiPhysics' `UniformUnsteady` active stress
exactly: S_a = T(t) e_F ⊗ e_F, with T(t) built from the same samples by the
same Fourier construction (`FourierTension`, ported from
FourierInterpolation.cpp; checked against an independent Python port in
`test_prescribed_active_tension_matches_closed_form`).

`coupled_3d.py` writes the 3D case in which the eighth-model cylinder replaces
the block, coupled to svZeroDSolver (`svZeroDSolver_interface`, `Coupled` Neu
BC on the endo face):

- The block is replaced by a cavity junction (zero-resistance stubs to the
  valves) whose extra inlet is the external FLOW block of the endo face.
  svMultiPhysics passes Q = ∫ v·n dA with the normal out of the wall, i.e.
  Q = −dV/dt: flow out of the cavity.
- The eighth model carries 1/8 of the flow, so the coupled circulation has
  all resistances and inductances ×8 and compliances ÷8, which reproduces the
  full-model pressures exactly.
- Coupling `implicit` with the linear solver GMRES: fsils' BiCGStab leaves out
  the coupled-BC resistance term of the tangent, which turns Newton into a
  fixed-point iteration that diverges once the circulation is stiffer than the
  wall (closed valves); `semi-implicit` converged only linearly. The
  finite-difference step of dP/dQ is set to SI flow scales (absolute 1e-14).
- The cavity pressure and the valve flows come from `results/svZeroD_data`,
  the 0D state that svMultiPhysics writes every step.
- The Linux interface library is built from this repository inside the Docker
  image into `build-docker/` (command in `coupled_3d.py`) and mounted at
  `/svzerod`.

`run_active_study.py` compares a cycle (volume, pressure, through-wall
displacement, stroke volume) with the Richardson-extrapolated 0D block in
the same circulation. The 0D model is periodic from the second cycle; the
first cycle, which starts from the same state in both models, is an equally
valid comparison without its first 0.02 s, in which the initial state
(reference volume at the initial pressure) is equilibrated
(`--cycles 1 --compare-cycle 1` halves the cost). Linear solver: fsils
GMRES up to nr = 2; Trilinos GMRES with algebraic multigrid
(`--ls trilinos-ml`, ~180 iterations per solve) from nr = 4, where fsils'
diagonal preconditioner stalls. Four processes are fastest on a 10-core
laptop; eight are slower (weaker multigrid, oversubscription).

`active_helical_L2/` is `active_prescribed/` with the helix of
`tests/cases/chamber_cylinder.json` (+60° to −60°) in the 2 L0 tube (initial
volume doubled, the same circulation; ejection fraction 0.18, peak 157 mmHg).
As in `passive_helical_L2/` the 3D model is the full cylinder over the full
length with free ends (u_z = 0, end dashpots), without circulation scaling,
and it needs `--ls trilinos-ml` at every level. Volume, pressure and flows are
those of the whole tube (they drive the circulation); the wall kinematics,
strains and the twist rate β are those of the mid-length section. The
contraction twists the 0D tube by up to 15.6° over its length, about 90× the
passive twist.

## Results

Maximum relative errors vs. the 0D reference; "order" is log2 of the error
ratio between successive meshes. Run times with 4 processes.

Passive inflation, 0–120 mmHg in 24 steps:

| Case | nr | V | ρ(R) | wall-thickness change | η | non-homogeneity | run time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `passive_isotropic` | 1 | 5.1e-3 | 2.1e-2 | 3.0e-2 | 3.2e-3 | 4.8e-2 | 3 s |
| | 2 | 4.1e-4 (3.6) | 2.8e-3 (2.9) | 6.6e-3 | 1.2e-3 | 1.1e-2 | 9 s |
| | 4 | 3.7e-5 (3.5) | 4.6e-4 (2.6) | 1.3e-3 | 2.5e-4 | 2.1e-3 | 101 s |
| `passive_anisotropic` | 1 | 2.6e-3 | 1.9e-2 | 3.0e-2 | 3.9e-3 | 4.7e-2 | 4 s |
| | 2 | 2.3e-4 (3.5) | 2.9e-3 (2.7) | 6.9e-3 | 1.4e-3 | 1.3e-2 | 12 s |
| | 4 | 2.1e-5 (3.5) | 4.9e-4 (2.6) | 1.4e-3 | 2.9e-4 | 2.5e-3 | 69 s |
| `passive_anisotropic_full` | 2 | 1.7e-4 | 1.6e-3 | 4.7e-3 | 1.2e-3 | 8.2e-3 | 36 s |
| | 4 | 1.5e-5 (3.5) | 2.8e-4 (2.6) | 9.0e-4 | 2.3e-4 | 1.6e-3 | 105 s |
| | 8 | 1.0e-6 (3.9) | 4.1e-5 (2.8) | 1.4e-4 | 3.3e-5 | 2.5e-4 | 40 min |
| | 16 | 7.6e-8 (3.8) | 5.5e-6 (2.9) | 1.8e-5 | 3.8e-6 | 4.0e-5 | 10.3 h (cluster, ML np 2) |

nr = 16 (149k elements) is beyond the fsils solvers on a laptop; it ran on the
cluster with Trilinos ML (see Slurm cluster). Its volume difference is the
same against the best estimate of the converged 0D solution (Richardson
128/256: 7.62e-8), so the 32/64 reference (within 7.8e-9 of it) does not
limit it; the direct 0D solution with ne = 256 (2.7e-7 off) would: against
it the nr = 16 volume difference reads 2.3e-7 (`output/summary_ne256.json`).

Strains and stresses through the wall over the same ramp (see Compared
quantities; J and σ_rr as element means, σ_rr relative to P_max):

| Case | nr | E_RR | E_TT | J | λ_f | σ_vm | σ_tt − σ_rr | σ_rr | σ_rr(Ri) + P |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `passive_isotropic` | 1 | 4.9e-2 | 3.6e-2 | 5.1e-1 | 3.4e-2 | 5.6e-2 | 5.7e-2 | 3.6e-1 | 3.8 |
| | 2 | 1.5e-2 (1.7) | 7.7e-3 (2.2) | 2.0e-1 (1.4) | 7.3e-3 (2.2) | 1.6e-2 (1.8) | 1.7e-2 | 1.3e-1 (1.4) | 1.4 |
| | 4 | 5.1e-3 (1.6) | 2.6e-3 (1.6) | 6.2e-2 (1.7) | 2.4e-3 (1.6) | 5.4e-3 (1.6) | 5.6e-3 | 4.3e-2 (1.7) | 4.6e-1 |
| `passive_anisotropic` | 1 | 6.1e-2 | 3.1e-2 | 4.6e-1 | 3.0e-2 | 5.7e-2 | 5.9e-2 | 2.4e-1 | 2.6 |
| | 2 | 2.0e-2 (1.6) | 7.2e-3 (2.1) | 1.6e-1 (1.5) | 7.0e-3 (2.1) | 1.7e-2 (1.7) | 1.8e-2 | 8.2e-2 (1.6) | 9.9e-1 |
| | 4 | 6.8e-3 (1.6) | 2.4e-3 (1.6) | 4.8e-2 (1.7) | 2.3e-3 (1.6) | 5.6e-3 (1.6) | 5.9e-3 | 2.4e-2 (1.7) | 3.2e-1 |
| `passive_anisotropic_full` | 2 | 1.4e-2 | 5.2e-3 | 3.0e-2 | 5.0e-3 | 3.3e-2 | 3.2e-2 | 3.4e-2 | 5.6e-1 |
| | 4 | 4.5e-3 (1.6) | 1.6e-3 (1.7) | 7.7e-3 (2.0) | 1.5e-3 (1.7) | 1.0e-2 (1.7) | 9.9e-3 | 8.8e-3 (1.9) | 1.7e-1 |
| | 8 | 1.3e-3 (1.8) | 4.7e-4 (1.8) | 1.8e-3 (2.1) | 4.4e-4 (1.8) | 2.7e-3 (1.9) | 2.7e-3 | 1.9e-3 (2.2) | 4.8e-2 |
| | 16 | 3.6e-4 (1.9) | 1.3e-4 (1.9) | 7.3e-4 (1.3) | 1.2e-4 (1.9) | 7.2e-4 (1.9) | 7.1e-4 | 7.0e-4 (1.5) | 1.3e-2 |

All converge, at order 1.6–2.2: they are gradients of the displacement,
recovered at the nodes by averaging the element gradients, and the largest
error sits at the endocardium. E_ZZ stays below 2e-3 of the in-plane
strains (the ends keep u_z = 0). The nodal hydrostatic stress oscillates
between the vertex and mid-edge levels of an element (κ (J − 1) of the nodal
J); its element means converge, but the nodal traction check σ_rr(Ri) = −P
is off by 0.3–0.5 P_max at nr = 4 with κ = 100 μ, and by 0.05 P_max at
nr = 8 and 0.013 P_max at nr = 16 in the stiffer full Holzapfel-Ogden wall. The 0D block, evaluated the
same way (nodal profiles, spline derivatives), converges at second order or
faster to the same reference (ne = 32: strains and deviatoric stresses
2e-5 to 3e-4, σ_rr(Ri) + P ≤ 3e-3).

Wall thickness and large deformation (eighth model): `passive_full_thin`
(Re/Ri = 1.2) and `passive_full_thick` (Re/Ri = 2.5) are the full
Holzapfel-Ogden case with the same inner radius, length and material;
`passive_isotropic_soft` is `passive_isotropic` with half the shear modulus
(μ = 50 kPa): the volume grows by 96% at 120 mmHg, 58% of the limit pressure
μ ln(Re/Ri) of the incompressible tube. The mesh keeps nr elements through
the wall and about cubic elements.

| Case | nr | elements | V | ρ(R) | wall-thickness change | E_RR | σ_vm | run time |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `passive_full_thin` | 1 | 432 | 2.0e-5 | 1.9e-4 | 1.9e-3 | 3.6e-3 | 6.7e-3 | 1 min |
| | 2 | 3456 | 1.4e-6 (3.8) | 3.3e-5 (2.5) | 4.0e-4 | 1.2e-3 | 1.6e-3 | 11 min |
| | 4 | 26040 | 1.5e-7 (3.2) | 5.2e-6 (2.7) | 6.9e-5 | 3.4e-4 | 3.8e-4 | 2.0 h (np 2) |
| `passive_full_thick` | 1 | 24 | 1.5e-2 | 7.9e-2 | 9.2e-2 | 1.6e-1 | 2.7e-1 | 3 s |
| | 2 | 144 | 1.8e-3 (3.1) | 1.2e-2 (2.7) | 1.8e-2 | 4.9e-2 | 8.7e-2 | 28 s |
| | 4 | 960 | 1.8e-4 (3.3) | 2.1e-3 (2.6) | 3.6e-3 | 1.7e-2 | 3.1e-2 | 5 min |
| | 8 | 6480 | 1.6e-5 (3.5) | 3.3e-4 (2.7) | 6.0e-4 | 5.4e-3 | 1.0e-2 | 20 min |
| | 16 | 48960 | 1.1e-6 (3.8) | 4.7e-5 (2.8) | 8.6e-5 | 1.5e-3 | 2.9e-3 | 3.2 h |
| `passive_isotropic_soft` | 1 | 54 | 5.3e-3 | 8.8e-3 | 1.4e-2 | 1.8e-2 | 3.9e-2 | 5 s |
| | 2 | 360 | 4.1e-4 (3.7) | 1.2e-3 (2.9) | 3.5e-3 | 5.6e-3 | 1.2e-2 | 27 s |
| | 4 | 2592 | 3.7e-5 (3.5) | 2.0e-4 (2.5) | 7.4e-4 | 1.9e-3 | 4.0e-3 | 5 min (np 1) |

Both walls converge at the rates of the base geometry. At the same nr the
differences grow with the wall thickness (larger elements, steeper
through-wall gradients): volume at nr = 4 1.5e-7 (thin), 1.5e-5 (base),
1.8e-4 (thick). The large deformation changes nothing: the soft wall's
volume errors are within 3% of those of `passive_isotropic`, its ρ errors
2.4× smaller relative to the larger displacement. Softer still, the
0D block stops near the limit pressure (μ = 30 and 24 kPa: at 93% and 96%
of μ ln(Re/Ri)), where no equilibrium exists under pressure control. The
smallest meshes need one process with fsils (BICG breaks down on 24
elements split four ways); thin nr = 4, thick nr = 16 and soft nr = 8 stall
in fsils BICG in the first step. The thin and thick tables are runs with
`--ls trilinos-ml` (thin np 2, thick np 4; the run times are with other
runs sharing the laptop, and multigrid is slow on the smallest meshes: thick
nr = 8 takes 3 min with fsils, 20 min with ML). Soft nr = 8 with ML needs
~8 min per step and was stopped.

Free length (closed ends, 6 L0 tube, eighth model with an end plate),
compared at the mid-length plane (V = π (Ri + ρ(Ri))² (1 + ε) L as in the
0D); nr = 4 and 8 ran on the cluster (Trilinos ML, np 2):

| Case | nr | elements | V | ρ(R) | wall-thickness change | ε | E_RR | E_TT | σ_vm | run time |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `passive_isotropic_closed` | 1 | 279 | 9.2e-03 | 3.4e-02 | 5.1e-02 | 8.1e-02 | 6.4e-02 | 5.2e-02 | 7.3e-02 | 14 s |
|  | 2 | 2124 | 1.3e-03 (2.8) | 4.6e-03 (2.9) | 1.0e-02 | 2.0e-02 (2.0) | 1.9e-02 | 1.0e-02 | 2.0e-02 | 6 min |
|  | 4 | 17136 | 2.1e-04 (2.6) | 7.5e-04 (2.6) | 2.0e-03 | 3.9e-03 (2.4) | 6.4e-03 | 3.1e-03 | 6.5e-03 | 38 min |
| `passive_anisotropic_closed` | 1 | 279 | 5.7e-03 | 3.3e-02 | 4.9e-02 | 3.6e-02 | 8.5e-02 | 4.9e-02 | 7.8e-02 | 23 s |
|  | 2 | 2124 | 9.1e-04 (2.6) | 5.2e-03 (2.7) | 1.0e-02 | 6.3e-03 (2.5) | 2.7e-02 | 1.0e-02 | 2.2e-02 | 6 min |
|  | 4 | 17136 | 1.6e-04 (2.5) | 9.1e-04 (2.5) | 1.9e-03 | 1.0e-03 (2.6) | 8.8e-03 | 3.1e-03 | 7.1e-03 | 31 min |
| `passive_anisotropic_full_closed` | 2 | 2124 | 7.7e-04 | 3.0e-03 | 6.1e-03 | 2.9e-03 | 1.9e-02 | 7.4e-03 | 4.2e-02 | 3 min |
|  | 4 | 17136 | 1.3e-04 (2.5) | 5.3e-04 (2.5) | 1.1e-03 | 3.7e-04 (3.0) | 5.9e-03 | 2.1e-03 | 1.2e-02 | 34 min |
|  | 8 | 135360 | 2.0e-05 (2.7) | 8.6e-05 (2.6) | 1.7e-04 | 3.0e-05 (3.6) | 1.7e-03 | 5.8e-04 | 3.2e-03 | 3.5 h |

All quantities converge to the 0D block, ρ at orders 2.5–2.9 and ε at 2.0–3.6;
the full Holzapfel-Ogden wall reaches 8.6e-5 in ρ and 3.0e-5 in ε at nr = 8.
The ρ differences are 1.6–2.1× those of the fixed-length tubes at the same
nr (the comparison is at the symmetry plane, see below). The volume follows
from ρ(Ri) here, not from the lumen surface, so it converges at the order of
ρ instead of 3.5. Inside the tube, 8–10 H from the plate, the differences are
smaller than at the mid-length plane (isotropic, nr = 4: endocardial ρ 1.0e-3
against 1.2e-3, ε 2–7e-4 against 3.9e-3): the symmetry plane, where every
node keeps u_z = 0, carries a larger discretization error, so the comparison
there is conservative. Strains and stresses converge at 1.6–2.3.

Tube length (isotropic wall): the 2 L0 tube puts the mid-length plane 4.1 H
from the plate, inside the end layer, and does not converge (ε 15–17% off);
6 L0 and 8 L0 agree to three digits. The 4 L0 tube at nr = 4 was stopped at
the 30 min limit of the local run.

| Tube | nr | elements | V | ρ(R) | wall-thickness change | ε | run time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 L0 | 1 | 135 | 6.9e-03 | 2.2e-02 | 4.7e-02 | 2.5e-01 | 7 s |
|  | 2 | 972 | 1.5e-03 | 1.3e-02 | 3.7e-03 | 1.7e-01 | 13 s |
|  | 4 | 7632 | 2.6e-03 | 1.4e-02 | 5.4e-03 | 1.5e-01 | 5 min |
| 4 L0 | 1 | 207 | 9.1e-03 | 3.4e-02 | 5.1e-02 | 7.5e-02 | 6 s |
|  | 2 | 1548 | 1.3e-03 | 4.7e-03 | 1.0e-02 | 1.8e-02 | 52 s |
| 6 L0 | 1 | 279 | 9.2e-03 | 3.4e-02 | 5.1e-02 | 8.1e-02 | 14 s |
|  | 2 | 2124 | 1.3e-03 | 4.6e-03 | 1.0e-02 | 2.0e-02 | 6 min |
|  | 4 | 17136 | 2.1e-04 | 7.5e-04 | 2.0e-03 | 3.9e-03 | 38 min |
| 8 L0 | 1 | 351 | 9.2e-03 | 3.4e-02 | 5.1e-02 | 8.1e-02 | 36 s |
|  | 2 | 2700 | 1.3e-03 | 4.6e-03 | 1.0e-02 | 2.0e-02 | 9 min |

Helical fibers (2 L0 tube, full cylinder with free ends), compared at the
mid-length section only (the middle element layer; β is the local twist rate
there, V = π (Ri + ρ(Ri))² L as in the 0D), max relative differences over the
ramp:

| Case | nr | elements | V | ρ(R) | wall thickness | β | φ | E_RR | E_TT | σ_vm | run time |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `passive_helical_L2` (±60°) | 1 | 648 | 1.5e-2 | 8.0e-2 | 9.6e-2 | 0.15 | 0.95 | 1.1e-1 | 1.0e-1 | 1.2e-1 | 3 min |
| | 2 | 4896 | 1.1e-3 (3.7) | 5.9e-3 (3.8) | 1.3e-2 (2.8) | 8.4e-2 (0.8) | 0.15 (2.6) | 2.6e-2 (2.0) | 1.4e-2 (2.8) | 2.7e-2 (2.2) | 12 min (np 8) |
| | 4 | 38016 | 1.7e-4 (2.8) | 8.7e-4 (2.8) | 3.1e-3 (2.1) | 1.6e-2 (2.4) | 1.6e-2 (3.3) | 8.2e-3 (1.6) | 3.9e-3 (1.9) | 8.3e-3 (1.7) | 5.0 h (np 2) |
| `passive_helical_a30_L2` (30°) | 1 | 648 | 6.7e-3 | 3.3e-2 | 4.9e-2 | 3.8e-2 | 1.8e-2 | 7.0e-2 | 5.6e-2 | 8.5e-2 | 2 min |
| | 2 | 4896 | 9.0e-4 (2.9) | 4.5e-3 (2.9) | 1.0e-2 (2.3) | 2.9e-3 (3.7) | 3.3e-3 (2.4) | 2.3e-2 (1.6) | 1.3e-2 (2.1) | 2.5e-2 (1.8) | 1.8 h (local) |
| `passive_helical_a60to0_L2` (+60° → 0°) | 1 | 648 | 7.0e-3 | 3.3e-2 | 5.0e-2 | 3.1e-3 | 2.4e-2 | 7.0e-2 | 5.6e-2 | 7.9e-2 | 6 min (local) |
| | 2 | 4896 | 8.9e-4 (3.0) | 4.2e-3 (3.0) | 1.1e-2 (2.2) | 1.9e-4 (4.0) | 4.3e-3 (2.5) | 2.3e-2 (1.6) | 1.3e-2 (2.1) | 2.4e-2 (1.7) | 50 min (local) |

(Trilinos ML at every level; nr = 1 and 4 of the ±60° helix and nr = 1 of the
30° helix on the cluster with 2 ranks; the same nr = 4 run on 1 rank gives the
same displacements to 2e-8.) Every quantity converges to the 0D at the middle
of the tube. Signed at 120 mmHg (3D/0D − 1 of the endocardial ρ; twist rate
3D/0D), the ±60° helix goes from −8.0e-2 / 1.110 (nr = 1) over −5.9e-3 / 1.084
to −8.7e-4 / 1.016 (nr = 4). Its twist rate converges late in relative terms:
the twist (−0.17° over 2 L0 at 120 mmHg) is a small difference of the
torques of the two fiber layers, which amplifies the stress error. The uniform
30° helix (−6.9° over 2 L0) and the helix from +60° to 0° (−5.6°) have large
twists, and there the twist rate converges at about fourth order: the torsion
of the 0D block is right, and the slow convergence of the ±60° twist rate is
the cancellation, not a model difference.

Active contraction (`active_prescribed`), first cycle, 400 steps per cycle
(P relative to the 158 mmHg range):

| nr | V | P | ρ(R) | stroke volume | run time |
| --- | --- | --- | --- | --- | --- |
| 1 | 1.3e-2 | 9.3e-3 | 4.9e-2 | 3.6e-2 | 2 min |
| 2 | 1.8e-3 (2.9) | 2.2e-3 (2.1) | 8.7e-3 (2.5) | 4.3e-3 (3.1) | 8 min |
| 4 | 3.9e-4 (2.2) | 2.3e-3 | 2.0e-3 (2.1) | 4.4e-4 (3.3) | 24 min |
| 8 | 3.9e-4 | 2.3e-3 | 1.5e-3 | 4.3e-5 (3.3) | 3.2 h (cluster, ML np 4) |

With half the time step (`active_prescribed_dt2`, 800 steps per cycle in both
models; run on the cluster):

| nr | V | P | ρ(R) | stroke volume | Q_in | valve times | run time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.3e-2 | 9.1e-3 | 4.9e-2 | 3.6e-2 | 5.3e-2 | 2.2 ms | 3 min |
| 2 | 1.6e-3 (3.0) | 1.2e-3 (3.0) | 8.2e-3 (2.6) | 4.3e-3 (3.1) | 6.8e-3 | 0.6 ms | 7 min |
| 4 | 2.2e-4 (2.9) | 8.8e-4 | 1.5e-3 (2.5) | 4.3e-4 (3.3) | 3.7e-3 | 0.8 ms | 50 min |

svMultiPhysics writes each row of `results/svZeroD_data` with the start time
of its time step but the 0D state at its end (`write_svZeroD_solution` is
called with `svZeroDTime` before that is advanced by dt, in
`svZeroD_subroutines.cpp`); `coupled_3d.read_svzerod_data` adds dt. Without
that correction the pressure error looked like a one-step lag in isovolumic
relaxation (−2300 mmHg/s): 3.6e-2 at nr = 4, halving with the time step (the
P column of the half-step runs above had the same offset, and those runs
were not kept). With it, the pressure error falls to 2.2e-3 at nr = 2 and
then stalls with the valve flows and dP/dt at the O(dt) difference between
the coupled 3D (generalized-α) and the stand-alone 0D time integration: half
the time step brings the nr = 4 pressure error from 2.3e-3 to 8.8e-4 and the
inflow error from 1.4e-2 to 3.7e-3, and the valve times then agree to less
than one step. The second cycle (`active_prescribed_c2`, `--cycles 2
--compare-cycle 2`), in which the 0D model is periodic, gives the errors of
the first (nr = 2: V 1.8e-3, P 1.3e-3, stroke volume 4.3e-3). The volume error
also falls with the time step around end-systole and early filling: with
half the step, volume, wall displacement and stroke volume converge at
third, 2.5th and third order. In the quasi-static late diastole the volume
error falls at fourth order (nr = 2 → 4: 5.5e-5 → 3.4e-6), as in the passive
cases.

Further quantities of the first cycle, 400 steps per cycle (flows relative to
their peak, valve times as the largest offset of the four openings and
closings):

| nr | Q_in | Q_out | stroke work | EF | max dP/dt | min dP/dt | valve times |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 5.2e-2 | 4.0e-2 | 4.4e-2 | 3.5e-2 | 2.9e-3 | 8.8e-3 | 2.4 ms |
| 2 | 1.4e-2 | 6.2e-3 | 5.7e-3 (2.9) | 4.2e-3 (3.1) | 1.0e-3 | 2.4e-3 | 1.5 ms |
| 4 | 1.4e-2 | 5.8e-3 | 1.1e-3 (2.4) | 4.3e-4 (3.3) | 8.6e-4 | 1.5e-3 | 1.6 ms |
| 8 | 1.4e-2 | 5.7e-3 | 5.8e-4 (0.9) | 4.4e-5 (3.3) | 7.9e-4 | 1.4e-3 | 1.6 ms |

| nr | E_RR | E_TT | J | λ_f |
| --- | --- | --- | --- | --- |
| 1 | 1.6e-1 | 4.7e-2 | 1.7e-1 | 6.4e-2 |
| 2 | 5.0e-2 (1.6) | 9.4e-3 (2.3) | 6.1e-2 (1.5) | 1.3e-2 (2.3) |
| 4 | 1.7e-2 (1.5) | 3.1e-3 (1.6) | 1.8e-2 (1.8) | 4.2e-3 (1.6) |
| 8 | 5.7e-3 (1.6) | 1.3e-3 (1.2) | 4.2e-3 (2.1) | 1.8e-3 (1.2) |

Stroke volume, ejection fraction and the wall strains converge with the mesh.
The valve flows, dP/dt and the valve times stall from nr = 2 on, like the
pressure; the valve times differ by less than one step (2.5 ms). At nr = 8
(cluster) the maxima of the volume, pressure and ρ differences stay at their
nr = 4 values. Split by cycle windows (as for the thick wall below):

| quantity, run | 0.02–0.1 s | 0.1–0.2 s | 0.2–0.35 s | 0.35–0.5 s | 0.5–1 s |
| --- | --- | --- | --- | --- | --- |
| V, nr 4 | 3.9e-4 | 2.1e-4 | 1.6e-4 | 3.7e-4 | 3.0e-5 |
| V, nr 8 | 3.9e-4 | 2.0e-4 | 6.3e-5 | 2.3e-4 | 3.1e-5 |
| V, nr 4, dt/2 | 1.9e-4 | 7.2e-5 | 1.5e-4 | 2.2e-4 | 6.6e-6 |
| ρ, nr 4 | 1.5e-3 | 8.0e-4 | 1.2e-3 | 1.9e-3 | 1.2e-4 |
| ρ, nr 8 | 1.5e-3 | 7.8e-4 | 2.5e-4 | 8.9e-4 | 1.2e-4 |
| ρ, nr 4, dt/2 | 7.3e-4 | 2.9e-4 | 1.2e-3 | 1.4e-3 | 6.2e-5 |

In ejection (0.2–0.35 s) the volume and ρ differences fall with the mesh and
not with the time step, in relaxation (0.35–0.5 s) with both, and early in
the cycle (0.02–0.2 s) with the time step and not with the mesh. At nr = 8 the maximum of all
three sits at the start of the compared cycle (t = 0.02 s): the remaining
difference is the time integration. Stroke work contains the pressure and
stalls with it (order 0.9 from nr = 4 to 8).

Thick wall (`active_thick`: Re/Ri = 2.5 instead of 1.74, the same inner
radius, material, tension and circulation; EF 0.37 instead of 0.27, peak
155 mmHg), first cycle, 400 steps per cycle:

| nr | elements | V | P | ρ(R) | stroke volume | stroke work | EF | Q_in | valve times | run time |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 24 | 1.2e-1 | 4.8e-2 | 2.5e-1 | 2.1e-1 | 2.4e-1 | 2.0e-1 | 2.3e-1 | 9.9 ms | 38 s (np 1) |
| 2 | 144 | 2.9e-2 (2.0) | 1.4e-2 (1.8) | 6.9e-2 (1.9) | 4.9e-2 (2.1) | 5.8e-2 (2.1) | 4.9e-2 (2.1) | 6.5e-2 | 1.6 ms | 5 min (np 1) |
| 4 | 960 | 4.5e-3 (2.7) | 2.4e-3 (2.5) | 1.3e-2 (2.4) | 7.0e-3 (2.8) | 8.6e-3 (2.8) | 7.0e-3 (2.8) | 1.8e-2 | 1.5 ms | 51 min (np 1) |
| 8 | 6480 | 7.2e-4 (2.6) | 2.3e-3 | 2.4e-3 (2.4) | 7.2e-4 (3.3) | 1.2e-3 (2.9) | 7.2e-4 (3.3) | 1.8e-2 | 1.6 ms | 2.4 h (ML, np 2) |

With half the time step (`active_thick_dt2`, 800 steps per cycle):

| nr | V | P | ρ(R) | stroke volume | Q_in | max dP/dt | valve times | run time |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.2e-1 | 4.8e-2 | 2.5e-1 | 2.1e-1 | 2.3e-1 | 1.6e-2 | 9.8 ms | 18 s (np 1) |
| 2 | 2.9e-2 | 1.4e-2 | 6.8e-2 | 4.9e-2 | 6.5e-2 | 2.7e-3 | 1.5 ms | 4 min (np 1) |
| 4 | 4.3e-3 | 2.3e-3 | 1.2e-2 | 7.0e-3 | 1.1e-2 | 1.0e-3 | 0.7 ms | 31 min (np 2) |

Volume, stroke volume, stroke work, ejection fraction and wall displacement
converge as in the base geometry. The pressure difference is about 2.3e-3 at
nr = 4 and 8, but for different reasons (maximum over sub-intervals of the
cycle, relative to the pressure range):

| run | 0.02–0.1 s | 0.1–0.2 s | 0.2–0.35 s | 0.35–0.5 s | 0.5–1 s |
| --- | --- | --- | --- | --- | --- |
| base nr 4 | 2.3e-3 | 9.5e-4 | 2.8e-4 | 1.1e-3 | 1.5e-4 |
| base nr 4, dt/2 | 8.8e-4 | 3.3e-4 | 1.5e-4 | 3.0e-4 | 2.3e-5 |
| thick nr 4 | 2.3e-3 | 1.5e-3 | 2.4e-3 | 2.3e-3 | 2.0e-4 |
| thick nr 4, dt/2 | 8.6e-4 | 6.7e-4 | 2.3e-3 | 2.2e-3 | 1.1e-4 |
| thick nr 8 | 2.3e-3 | 1.4e-3 | 5.0e-4 | 1.7e-3 | 2.2e-4 |

At the start of the compared cycle (0.02 s) the difference is the time
integration in both walls (halved time step: 2.3e-3 → 8.7e-4); in the thick
wall at nr = 4 a second, spatial part of the same size sits in the ejection
and relaxation (unchanged with half the step), which falls at nr = 8. Half the
step also reduces the inflow difference (1.8e-2 → 1.1e-2), the max dP/dt
difference (2.4e-3 → 1.0e-3) and the valve-time offset (1.5 → 0.7 ms).
(nr = 8 shared the laptop with other runs.)

Helical fibers (`active_helical_L2`, 2 L0 tube), first cycle, 400 steps per
cycle (ρ, strains and β at mid-length; β relative to its peak):

| nr | V | P | ρ(R) | stroke volume | stroke work | EF | β | Q_in | Q_out | valve times | run time |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.7e-2 | 5.7e-2 | 1.1e-1 | 4.5e-2 | 6.2e-2 | 4.7e-2 | 4.6e-2 | 4.2e-1 | 2.3e-1 | 10.6 ms | 10 min (cluster, np 2) |
| 2 | 5.5e-3 | 1.2e-2 | 2.5e-2 | 2.6e-2 | 3.3e-2 | 2.6e-2 | 6.9e-3 | 1.1e-1 | 4.5e-2 | 0.9 ms | 2.0 h (cluster, np 2) |

From nr = 1 to 2 every difference falls. The stroke volume falls least: the
circulation sees the volume of the whole tube, end layers included. The 0D
block loaded with the 3D cavity pressure (`zero_d_replay.py`, see 0D
reference and summary) separates the two: at nr = 2 the stroke volume of the
whole tube differs from the 0D by 1.5e-2, that of the mid-length section
(V = π (Ri + ρ(Ri))² L) by 9e-4. In the test-case tube the end layers are
twice the share of the volume, and the pressure and stroke-volume differences
stayed at 3–6% from nr = 2 to 4; hence the 2 L0 tube.

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
| `eps` | axial strain ε (free length; 3D: fit of u_z at the mid-length plane) | relative to max \|ε_ref\| |

The reference is the closed form for the incompressible isotropic tube and
otherwise the Richardson extrapolation of the 0D block at `num_elements` and
twice that (the 0D model converges at second order). The study also reports
the 0D convergence against the same reference and the non-homogeneity of the
3D solution (largest deviation of ρ from its level mean).

Strains and stresses are compared as through-wall profiles: level means of the
3D nodal fields (below) vs the same quantities of the 0D solution
(`zero_d_fields`), interpolated to the 3D levels by cubic splines. The 0D
deformation gradient follows from the nodal ρ, φ, η and β, with the radial
derivatives from cubic splines of the nodal values (the 0D derivatives are
element-wise constant), and ε (0 with fixed length); its stresses come from the 0D material law
with the ST91 pressure of the local J.

| Key | Quantity | Error measure |
| --- | --- | --- |
| `E_RR`, `E_TT`, `E_ZZ` | Green-Lagrange strain, reference cylindrical basis | relative to max(\|E_RR\|, \|E_TT\|) |
| `J` | Jacobian, element means | relative to max \|J − 1\| |
| `lambda_f` | fiber stretch | relative to max \|λ_f − 1\| |
| `sigma_vm` | von Mises stress | relative to its maximum |
| `sigma_dev` | σ_tt − σ_rr (Cauchy, deformed basis) | relative to its maximum |
| `sigma_rr` | radial Cauchy stress, element means | relative to P_max |
| `sigma_rr_endo` | σ_rr(Ri) + P, the endo traction condition | relative to P_max |

Element means use Simpson's rule over the three node levels of a P2 element;
the 0D counterpart is the mean of its spline over the element.

`run_active_study.py` compares over the cycle the volume, pressure, ρ(R) and
stroke volume, the inlet and outlet valve flows (relative to their peak), the
valve opening and closing times (crossings of 2% of the peak flow, in ms),
the stroke work (PV-loop area), ejection fraction and the extrema of dP/dt,
and the strain, Jacobian and fiber-stretch profiles as above.

Each study writes `CASE/output/summary.json` and one multi-panel PDF:
`comparison_<formulation>.pdf` (passive: PV, ρ(R), strain, J, λ_f, stress
profiles at 120 mmHg and the convergence of all errors) or
`comparison_cycle<k>.pdf` (active: PV loop, volume, pressure, valve flows,
dP/dt, endo and epi fiber stretch, strain profiles at end systole and end
diastole, convergence).

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
smooth; J and σ_rr are therefore compared as element means. The coupled
active case adds the kinematic profiles (J, E_RR, E_TT, E_ZZ, λ_f) to
`coupled.json`.

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

Helical fibers need the full length; from nr = 2 on, the Trilinos solver:

```bash
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/run_study.py scripts/ChamberCylinder/verification_3d/passive_helical_L2 --nr 2 4 --run --full-length --ls trilinos-ml --np 2
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/helical_figure.py scripts/ChamberCylinder/verification_3d/passive_helical_L2
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/helical_sketch.py scripts/ChamberCylinder/verification_3d/passive_helical_L2
```

Free length (ends `lids` in the 0D formulation: closed by an end plate,
the eighth model, no extra option); the figure needs nr = 1, 2, 4 of the 6 L0
tube and nr = 2 of the 2, 4 and 8 L0 tubes:

```bash
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/run_study.py scripts/ChamberCylinder/verification_3d/passive_isotropic_closed --nr 1 2 4 --run --np 4
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/closed_end_figure.py
```

`--nr` sets the elements through the wall (the mesh-size parameter;
circumferential and axial divisions follow from h = (Re − Ri)/nr). Cases go to
`CASE/runs/`, results and the comparison PDF to `CASE/output/`, which also
caches the 0D runs. Without `--run` the script post-processes existing
results (as they were run) and generates the missing cases.

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

### 0D reference and summary

`zero_d_convergence.py` runs the 0D input of a case at several radial
resolutions (cached in `CASE/output/` like the studies) and reports how far
each level and each Richardson extrapolation is from the finest one
(`CASE/output/zero_d_convergence.json`). `--reference-ne N` of both studies
compares the 3D results with the direct 0D solution with N elements instead of
the Richardson extrapolation (`CASE/output/summary_ne<N>.json`, no figure).
The extrapolations converge at fourth order: the passive references (32/64)
are within 8e-9 of the extrapolation of the two finest levels in the volume
(3e-7 in the twist rate), while the direct ne = 256 solution is off by up to
2.7e-7 (3.7e-5) — enough to shift the nr = 8 volume differences by 5%, so the
extrapolated references are kept. `zero_d_replay.py` is a diagnostic for the
coupled cases: it loads the 0D block alone with the cavity pressure of a
coupled 3D run (`CASE/output/replay.json`), so that both are compared at the
same pressure, without the feedback of the circulation (the 3D model is
unchanged). The replayed 0D solution is the Richardson extrapolation of the
coupled reference (ne = 16/32; `--ne N` for a direct solution). `make_summary.py` collects all study summaries into
`summary/verification_summary.pdf` (pdflatex; the narrative is in
`summary/text/`).

```bash
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/zero_d_convergence.py scripts/ChamberCylinder/verification_3d/passive_anisotropic_full --ne 16 32 64 128 256
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/zero_d_convergence.py scripts/ChamberCylinder/verification_3d/active_prescribed --active --ne 16 32 64 128
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/run_study.py scripts/ChamberCylinder/verification_3d/passive_anisotropic_full --nr 2 4 8 16 --reference-ne 256
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/zero_d_replay.py scripts/ChamberCylinder/verification_3d/active_helical_L2 --nr 1 2
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/make_summary.py
```

### Slurm cluster

`cluster.py` runs generated cases on a Slurm cluster over SSH, in an
Apptainer image of the same Docker image (`cluster.py image` builds it with
`cluster/build_image.sbatch`). Log in once with SSH multiplexing
(`ControlMaster`/`ControlPersist` for the host), so that the script's
non-interactive SSH commands need no second-factor login, and set
`SVZEROD_CLUSTER` (SSH host), `SVZEROD_CLUSTER_WORK` (remote scratch folder)
and, if needed, `SVZEROD_CLUSTER_ACCOUNT`. Generate the cases locally (a study
without `--run`), then:

```bash
python scripts/ChamberCylinder/verification_3d/cluster.py submit CASE/runs/struct-p2_nr4 --np 4 --partition day --time 12:00:00
python scripts/ChamberCylinder/verification_3d/cluster.py watch CASE/runs/struct-p2_nr4
python scripts/ChamberCylinder/verification_3d/cluster.py fetch CASE/runs/struct-p2_nr4
```

and post-process locally as usual. `watch` streams the remote solver log
through the same checks as a local run and cancels a job on NaN or a
diverging residual; a time step that ends above the Newton tolerance (which
svMultiPhysics passes on to the next step) is reported, not cancelled. (The
strict local check once cancelled two 15-hour runs when a restarted watch
read their logs from the start.) The solver input keeps one restart file,
overwritten every step (`results/stFile_last.bin`), so an interrupted run
can be continued with `<Continue_previous_simulation> 1`. Coupled
cases need the Linux interface library in `build-docker/` (it is uploaded and
bound to `/svzerod`). On Yale's Bouchet the native runs are 5–6× faster than
the emulated local runs and give the same results (bitwise for coupled runs
with the same number of MPI ranks). Trilinos ML degrades with the number of
ranks (no residual reduction at 32 ranks, breakdowns at 8): use 4 ranks up to
about 20k elements and 1–2 for the large full cylinders. The image also has a
PETSc build of svMultiPhysics, but svMultiPhysics does not pass
`<Configuration_file>` on to PETSc, and options through `PETSC_OPTIONS`
crashed, so a direct solver (MUMPS) was not available.
