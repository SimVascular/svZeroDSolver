@page block_types Block Types in the Input File

[TOC]

Every block in an svZeroDSolver input file is created from a *type* string,
read from a key that depends on the section the block is in. This page lists
the valid type strings for each section and links to each block's
documentation, which describes its parameters and JSON usage.

If a type is misspelled, or a block is used in the wrong section, loading the
file fails with an error that names the block, the section, and the key to
fix, and lists the valid values. That list is generated from the blocks
registered in the solver, so it is the authoritative version of the tables
below.

# Type keys by section

| Section | Type key | Valid values |
|---|---|---|
| `vessels` | `zero_d_element_type` | `BloodVessel` (BloodVessel), `BloodVesselCRL` (BloodVesselCRL), `BloodVesselRC` (BloodVesselRC) |
| `junctions` | `junction_type` | `NORMAL_JUNCTION` (Junction), `resistive_junction` (ResistiveJunction), `BloodVesselJunction` (BloodVesselJunction) |
| `boundary_conditions` | `bc_type` | see [Boundary condition types](#bc_types) |
| `external_solver_coupling_blocks` | `type` | the same values as `bc_type`; in practice `FLOW` or `PRESSURE` |
| `closed_loop_blocks` | `closed_loop_type` | `ClosedLoopHeartAndPulmonary` (ClosedLoopHeartPulmonary) |
| `valves` | `type` | `ValveTanh` (ValveTanh), `PiecewiseValve` (PiecewiseValve) |
| `chambers` | `type` | `ChamberSphere` (ChamberSphere), `ChamberElastanceInductor` (ChamberElastanceInductor), `ChamberElastanceInductorExponential` (ChamberElastanceInductorExponential), `LinearElastanceChamber` (LinearElastanceChamber) |

Type strings are case-sensitive.

# Boundary condition types {#bc_types}

| `bc_type` | Block |
|---|---|
| `FLOW` | FlowReferenceBC: prescribed flow rate |
| `PRESSURE` | PressureReferenceBC: prescribed pressure |
| `RESISTANCE` | ResistanceBC |
| `RCR` | WindkesselBC |
| `IMPEDANCE` | ImpedanceBC |
| `CORONARY` | OpenLoopCoronaryBC |
| `CORONARY_VAR_RES` | OpenLoopCoronaryVarResBC |
| `ClosedLoopRCR` | ClosedLoopRCRBC |
| `ClosedLoopCoronaryLeft` | ClosedLoopCoronaryLeftBC |
| `ClosedLoopCoronaryRight` | ClosedLoopCoronaryRightBC |

`ClosedLoopHeartAndPulmonary` is not a `bc_type`: the heart is configured
under `closed_loop_blocks`.

# External coupling blocks

An entry in `external_solver_coupling_blocks` is the boundary through which an
external solver (for example a 3D solver such as svMultiPhysics) exchanges
flow and pressure with the 0D model. Besides `name`, `type`, and `values`, it
has:

- `connected_block`: the 0D block on the other side of the coupling block,
  given as a boundary condition's `bc_name`, a vessel's `vessel_name`, or
  `ClosedLoopHeartAndPulmonary`.
- `location`: `inlet` or `outlet`, the side of `connected_block` the coupling
  block attaches to. With `inlet`, flow goes from the coupling block into
  `connected_block`; with `outlet`, from `connected_block` into the coupling
  block.

Only some connected blocks are supported at each location. A connected
block's type is its `bc_type` for a boundary condition, or `BloodVessel` for
any vessel:

| `location` | Supported connected block types |
|---|---|
| `inlet` | `RESISTANCE`, `RCR`, `IMPEDANCE`, `ClosedLoopRCR`, `CORONARY`, `ClosedLoopCoronaryLeft`, `ClosedLoopCoronaryRight`, `BloodVessel` |
| `outlet` | `ClosedLoopRCR`, `ClosedLoopHeartAndPulmonary`, `BloodVessel` |
