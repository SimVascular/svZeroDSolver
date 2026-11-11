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
| 2c | `passive_helical/`: the wall of step 2 with helical fibers (+60° to −60°), full cylinder with free ends | implemented (nr 1, 2, 4) |
| 3a | `active_prescribed/`: PV loop with a prescribed active tension, 3D cylinder coupled to the circulation of `tests/cases/chamber_cylinder.json` | implemented (nr 1, 2) |
| 3a' | `active_helical/`: the same with helical fibers (+60° to −60°), full 3D cylinder with free ends | implemented (nr 1) |
| 3b | PV loop with an active-stress ODE model (e.g. Land-Niederer in both) | planned |

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

`passive_helical/` is the wall of `passive_anisotropic/` with the helix of
`tests/cases/chamber_cylinder.json`, +60° at the endocardium to −60° at the
epicardium. Mirroring θ or Z flips the helix angle, so the eighth model's
symmetry planes do not apply: the 3D model is the full cylinder over the full
length. It is as free as the 0D: both ends keep u_z = 0 (fixed length) and are
otherwise free to twist and expand. In-plane dashpots on both end faces
(Robin BCs with zero stiffness, damping 1e-3 μ/Re·dt) hold the three
rigid-body modes (x, y translation, rotation about z): within a load step they
act as springs of 1e-3 μ/Re, at rest they carry no force (during the ramp
< 2e-5 of the cavity pressure). No 0D data enters the 3D model. Alternatives
that did worse: a weak epicardial spring (1e-6 μ/Re) gives the same result but
175–950 linear iterations per solve; no support at all with GMRES stalls
(Newton stops at 4e-9, the discrete system is not exactly consistent); two
pinned radial lines at the base give the same averages but disturb the wall
near the pins (non-homogeneity 8.5% vs 3.4% at nr = 1).

The uniform 0D state does not solve this 3D problem exactly. Every 0D
cross-section carries the through-wall shear σ_θZ(R) (positive near the
endocardium, negative near the epicardium, zero net torque, up to 18% of
σ_θθ at 120 mmHg); the free 3D ends carry none. Near the ends the 3D wall
twists faster and expands up to 1.6% more. The tube is about four wall
thicknesses long (L = 57 mm, H = 14 mm), so for the small twist of the ±60°
helix (βL = −0.085° at 120 mmHg) the end layers reach the middle. This model
difference does not shrink with the mesh, only with a longer tube.
`run_study.py` therefore compares helical cases at the mid-length section
only (the middle element layer). `helical_sketch.py` draws both models and `helical_figure.py` the profiles
along the axis, with the 0D section shear and twist against the helix angle
(`output/helical_bcs_sketch.pdf`, `output/helical_discrepancy.pdf`).

The full cylinder is expensive: it has about 8× the elements of the eighth
model at the same nr, and it is held only weakly against rigid motion. The
fsils BiCGStab solver stalls from nr = 2 on, with any support; Trilinos GMRES
with ML (`--ls trilinos-ml`) converges with 500–1000 iterations per solve.
At nr = 2 (2592 elements, as many as the eighth model at nr = 4) the ramp takes
12 min, against 69 s for the eighth model; nr = 4 (19584 elements) takes
82 min with 4 MPI ranks.

### Active contraction in a circulation

`active_prescribed/` is the 0D test case `tests/cases/chamber_cylinder.json`
(pressure source, tanh valves, RC vessels, 1 s period) with the cylinder
adapted to what svMultiPhysics can represent: fixed length, circumferential
fibers, the Holzapfel-Ogden wall of `passive_anisotropic/`, no viscosity, no
wall inertia, and the BCS active law replaced by a prescribed tension
T(t) = 80 kPa sin²(π (t − 0.1 s)/0.4 s) on 0.1–0.5 s. The 0D block's
`active_tension` reproduces svMultiPhysics' `UniformUnsteady` active stress
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

`active_helical/` is `active_prescribed/` with the helix of
`tests/cases/chamber_cylinder.json` (+60° to −60°). As in `passive_helical/`
the 3D model is the full cylinder over the full length with free ends (u_z = 0,
end dashpots), without circulation scaling, and it needs `--ls trilinos-ml`
at every level. Volume, pressure and flows are those of the whole tube (they
drive the circulation); the wall kinematics, strains and the twist rate β are
those of the mid-length section. The contraction twists the 0D tube by up to
7.9° over its length, 100× the passive twist.

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
| | 8 | 1.1e-6 (3.9) | 4.1e-5 (2.8) | 1.4e-4 | 3.3e-5 | 2.5e-4 | 40 min |

nr = 16 (149k elements) is beyond the fsils solvers on a laptop.

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

All converge, at order 1.6–2.2: they are gradients of the displacement,
recovered at the nodes by averaging the element gradients, and the largest
error sits at the endocardium. E_ZZ stays below 2e-3 of the in-plane
strains (the ends keep u_z = 0). The nodal hydrostatic stress oscillates
between the vertex and mid-edge levels of an element (κ (J − 1) of the nodal
J); its element means converge, but the nodal traction check σ_rr(Ri) = −P
is off by 0.3–0.5 P_max at nr = 4 with κ = 100 μ, and by 0.05 P_max at
nr = 8 in the stiffer full Holzapfel-Ogden wall. The 0D block, evaluated the
same way (nodal profiles, spline derivatives), converges at second order or
faster to the same reference (ne = 32: strains and deviatoric stresses
2e-5 to 3e-4, σ_rr(Ri) + P ≤ 3e-3).

Helical fibers (`passive_helical`, full cylinder with free ends), compared at
the mid-length section only (the middle element layer; β is the local twist
rate there, V = π (Ri + ρ(Ri))² L as in the 0D), max relative differences
over the ramp:

| nr | V | ρ(R) | wall thickness | β | φ | E_RR | E_TT | E_ZZ | σ_vm | non-homogeneity | run time |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.4e-2 | 7.5e-2 | 8.9e-2 | 0.17 | 1.9 | 1.1e-1 | 9.7e-2 | 3.8e-3 | 1.2e-1 | 1.8e-2 | 2 min (fsils) |
| 2 | 4.9e-4 | 2.5e-3 | 1.2e-2 | 0.23 | 0.31 | 2.3e-2 | 1.1e-2 | 4.7e-3 | 2.4e-2 | 6.3e-3 | 12 min (Trilinos ML) |
| 4 | 4.8e-4 | 3.2e-3 | 1.7e-3 | 0.16 | 3.2e-2 | 5.5e-3 | 3.0e-3 | 4.1e-3 | 5.2e-3 | 1.4e-3 | 82 min (Trilinos ML) |

At mid-length the wall thickness, the transmural rotation, the strains and
stresses and the axisymmetry converge toward the 0D. The volume (5e-4), the
radial displacement (3e-3), the twist rate (16%) and the axial strain stall
from nr = 2 on: the tube is only about four wall thicknesses long, and the
end effect still reaches its middle. Along the whole tube
(`helical_discrepancy.pdf`, where the nr = 2 and 4 profiles hardly differ)
the endocardium expands 0.3% more than the 0D at mid-length and 1.6% more at
the ends, and the tube twists 1.9× as much over its length (−0.16° vs
−0.085° at 120 mmHg). With circumferential fibers the discretization error at
nr = 4 is about 1e-4, so these differences are the model difference between
the free-ended tube and the uniform 0D state, not mesh error.

Active contraction (`active_prescribed`), first cycle, 400 steps per cycle
(P relative to the 158 mmHg range):

| nr | V | P | ρ(R) | stroke volume | run time |
| --- | --- | --- | --- | --- | --- |
| 1 | 1.3e-2 | 9.3e-3 | 4.9e-2 | 3.6e-2 | 2 min |
| 2 | 1.8e-3 (2.9) | 2.2e-3 (2.1) | 8.7e-3 (2.5) | 4.3e-3 (3.1) | 8 min |
| 4 | 3.9e-4 (2.2) | 2.3e-3 | 2.0e-3 (2.1) | 4.4e-4 (3.3) | 24 min |

With half the time step (800 steps per cycle; P not listed, see below):

| nr | V | ρ(R) | stroke volume | run time |
| --- | --- | --- | --- | --- |
| 1 | 1.3e-2 | 4.9e-2 | 3.4e-2 | 4 min |
| 2 | 1.6e-3 (3.0) | 8.2e-3 (2.6) | 4.3e-3 (3.0) | 13 min |
| 4 | 2.2e-4 (2.9) | 1.5e-3 (2.5) | 4.3e-4 (3.3) | 44 min |

svMultiPhysics writes each row of `results/svZeroD_data` with the start time
of its time step but the 0D state at its end (`write_svZeroD_solution` is
called with `svZeroDTime` before that is advanced by dt, in
`svZeroD_subroutines.cpp`); `coupled_3d.read_svzerod_data` adds dt. Without
that correction the pressure error looked like a one-step lag in isovolumic
relaxation (−2300 mmHg/s): 3.6e-2 at nr = 4, halving with the time step (the
P column of the half-step runs above had the same offset, and those runs
were not kept). With it, the pressure error falls to 2.2e-3 at nr = 2 and
then stalls with the valve flows and dP/dt, presumably at the O(dt)
difference between the coupled 3D (generalized-α) and the stand-alone 0D
time integration (not yet confirmed by a half-step run). The volume error
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

| nr | E_RR | E_TT | J | λ_f |
| --- | --- | --- | --- | --- |
| 1 | 1.6e-1 | 4.7e-2 | 1.7e-1 | 6.4e-2 |
| 2 | 5.0e-2 (1.6) | 9.4e-3 (2.3) | 6.1e-2 (1.5) | 1.3e-2 (2.3) |
| 4 | 1.7e-2 (1.5) | 3.1e-3 (1.6) | 1.8e-2 (1.8) | 4.2e-3 (1.6) |

Stroke work, ejection fraction and the wall strains converge with the mesh.
The valve flows, dP/dt and the valve times stall from nr = 2 on, like the
pressure; the valve times differ by less than one step (2.5 ms).

Helical fibers (`active_helical`), first cycle, 400 steps per cycle (ρ, strains
and β at mid-length; β relative to its peak):

| nr | V | P | ρ(R) | stroke volume | stroke work | EF | β | Q_in | Q_out | valve times | run time |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.4e-2 | 4.7e-2 | 8.8e-2 | 1.6e-2 | 2.4e-2 | 1.8e-2 | 5.3e-2 | 3.2e-1 | 2.4e-1 | 8.0 ms | 8 min |

At nr = 1 the PV loop, pressure and dP/dt nearly coincide with the 0D, and the
mid-length twist follows the 0D's (peak 7.4° vs 7.9° over the length); the
flow errors are the steep valve-flow fronts shifted by up to 8 ms.

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

Strains and stresses are compared as through-wall profiles: level means of the
3D nodal fields (below) vs the same quantities of the 0D solution
(`zero_d_fields`), interpolated to the 3D levels by cubic splines. The 0D
deformation gradient follows from the nodal ρ, φ, η and β, with the radial
derivatives from cubic splines of the nodal values (the 0D derivatives are
element-wise constant), and ε = 0; its stresses come from the 0D material law
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
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/run_study.py scripts/ChamberCylinder/verification_3d/passive_helical --nr 2 4 --run --full-length --ls trilinos-ml --np 4
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/helical_figure.py scripts/ChamberCylinder/verification_3d/passive_helical
PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/helical_sketch.py scripts/ChamberCylinder/verification_3d/passive_helical
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
