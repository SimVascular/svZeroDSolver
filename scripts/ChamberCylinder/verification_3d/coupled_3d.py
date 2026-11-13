# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""3D svMultiPhysics cylinder coupled to the 0D circulation of a ChamberCylinder.

Reads an svZeroDSolver input in which a ChamberCylinder block sits in a
circulation (valves, vessels, pressure boundary conditions) and writes the
svMultiPhysics case in which the 3D cylinder replaces the block:

  CASE_DIR/solver.xml               struct, HolzapfelOgden (or neo-Hookean),
                                    UniformUnsteady active stress, endo face
                                    coupled to svZeroDSolver (semi-implicit)
  CASE_DIR/svzerod_3Dcoupling.json  the circulation, with the cylinder replaced
                                    by a cavity junction whose third branch is
                                    the external FLOW block of the endo face
  CASE_DIR/active_tension.dat       the tension samples of the 0D input
  CASE_DIR/mesh/                    eighth model (90-degree sector, half
                                    length), nodal fiber and sheet directions

The eighth model carries 1/8 of the cavity flow, so every resistance and
inductance of the coupled circulation is multiplied and every compliance
divided by 8: with the same pressures, each block then carries 1/8 of the
flow of the full model, which is exact for the linear vessels and for valves
that switch on the pressure difference.

The svZeroDSolver interface library (Linux, for the Docker image) is built
from this repository with
  docker run --rm --platform linux/amd64 -v $PWD:/svzerod -w /svzerod \\
      simvascular/solver:latest bash -c \\
      'cmake -S . -B build-docker -DCMAKE_BUILD_TYPE=Release && \\
       cmake --build build-docker --target svzero_interface -j4'
and mounted into the container at /svzerod (see run()).

Usage:
  python coupled_3d.py generate MODEL_0D.json CASE_DIR --nr 1 [--cycles 3]
  python coupled_3d.py run CASE_DIR [--np 4]
  python coupled_3d.py post CASE_DIR
"""
import argparse
import copy
import json
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cylinder_3d as C3  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
LIB_IN_IMAGE = "/svzerod/build-docker/src/interface/libsvzero_interface.so"
BLOCK = "LV"  # external coupling block of the endo face


def find_cylinder(cfg):
    blocks = [v for v in cfg["vessels"] if v["zero_d_element_type"] == "ChamberCylinder"]
    if len(blocks) != 1:
        raise ValueError(f"expected one ChamberCylinder vessel, found {len(blocks)}")
    return blocks[0]


def coupling_config(cfg, scale, n_substeps=2):
    """svZeroDSolver input of the circulation coupled to the 3D cavity."""
    c = copy.deepcopy(cfg)
    cyl = find_cylinder(c)
    name = cyl["vessel_name"]
    c["vessels"] = [v for v in c["vessels"] if v is not cyl]
    for v in c["vessels"]:
        vals = v["zero_d_element_values"]
        for k in ("R_poiseuille", "L"):
            if k in vals:
                vals[k] *= scale
        if "C" in vals:
            vals["C"] /= scale
        if "stenosis_coefficient" in vals:  # Delta P = K Q |Q|
            vals["stenosis_coefficient"] *= scale * scale
    for valve in c.get("valves", []):
        p = valve["params"]
        for k in ("Rmax", "Rmin"):
            p[k] *= scale
        if p["downstream_block"] == name:
            p["downstream_block"] = "cavity_in"
        if p["upstream_block"] == name:
            p["upstream_block"] = "cavity_out"
    vid = max(v["vessel_id"] for v in c["vessels"]) + 1
    for i, stub in enumerate(("cavity_in", "cavity_out", "cavity_3d")):
        c["vessels"].append({"vessel_id": vid + i, "vessel_length": 1.0, "vessel_name": stub,
                             "zero_d_element_type": "BloodVessel",
                             "zero_d_element_values": {"R_poiseuille": 0.0, "C": 0.0}})
    # svMultiPhysics passes the flux int v.n dA of the endo face with the normal
    # out of the wall, i.e. into the cavity: Q = -dV/dt, flow out of the cavity.
    # It enters the cavity junction as a second inlet: Q_in + Q = Q_out.
    c.setdefault("junctions", []).append(
        {"junction_name": "cavity", "junction_type": "NORMAL_JUNCTION",
         "inlet_blocks": ["cavity_in", "cavity_3d"], "outlet_blocks": ["cavity_out"]})
    c["external_solver_coupling_blocks"] = [
        {"name": BLOCK, "type": "FLOW", "location": "inlet", "connected_block": "cavity_3d",
         "periodic": False, "values": {"Q": [0.0, 0.0], "t": [0.0, 1.0]}}]
    sp = c["simulation_parameters"]
    for k in ("number_of_cardiac_cycles", "number_of_time_pts_per_cardiac_cycle",
              "output_interval", "output_variable_based"):
        sp.pop(k, None)
    sp.update(coupled_simulation=True, number_of_time_pts=n_substeps, output_all_cycles=True,
              steady_initial=False)
    ic = c.get("initial_condition", {})
    c["initial_condition"] = {k: v for k, v in ic.items() if not k.endswith(":" + name)}
    return c


def generate(model_0d, case_dir, nr, cycles=None, density=1.0, newton_tol=1e-10,
             linear_solver=None, n_substeps=2, end_damping=1e-3):
    """Write the coupled svMultiPhysics case of a ChamberCylinder circulation.

    Circumferential or axial fibers use the eighth model (the circulation
    impedances scaled by 8); helical fibers the full cylinder with in-plane
    dashpots on both ends (damping end_damping * mu / Re * dt, see
    cylinder_3d.generate).
    """
    with open(model_0d) as f:
        cfg = json.load(f)
    cyl = find_cylinder(cfg)
    v, form = C3.flat_values(cyl), C3.formulation(cyl)
    if form["ends"] != "fixed_length":
        raise ValueError('the 3D model has fixed-length ends: set the 0D formulation '
                         '"ends": "fixed_length"')
    at = cyl["active_stress"]
    if form["active"] != "uniform_unsteady":
        raise ValueError('the 3D model needs the 0D active_stress "uniform_unsteady" '
                         '(a prescribed tension, svMultiPhysics UniformUnsteady)')
    mat = C3.material(v, form)
    sp = cfg["simulation_parameters"]
    period = float(sp["cardiac_period"])
    n_pts = int(sp["number_of_time_pts_per_cardiac_cycle"])
    cycles = int(cycles or sp.get("number_of_cardiac_cycles", 1))
    # circumferential or axial fibers: the eighth model; helical fibers: the
    # full cylinder over the full length, with end dashpots (cylinder_3d)
    helical = C3.helical(mat)
    sector, Lm = (360.0, v["length"]) if helical else (90.0, 0.5 * v["length"])
    scale = (360.0 / sector) * v["length"] / Lm
    os.makedirs(case_dir, exist_ok=True)
    fm = C3.FORMULATIONS["struct-p2"]
    mesh = C3.mesh_cylinder(v["Ri"], v["Re"], Lm, nr, sector, fm["order"])
    C3.write_mesh(mesh, os.path.join(case_dir, "mesh"))
    from cylinder_fibers import cylinder_fibers, write_fibers
    fib = cylinder_fibers(mesh["points"], mesh["cells"], v["Ri"], v["Re"],
                          mat.get("alpha_endo", 0.0), mat.get("alpha_epi", 0.0), "points")
    paths = write_fibers(os.path.join(case_dir, "mesh"), mesh["points"], mesh["cells"], 24, fib)
    with open(os.path.join(case_dir, "active_tension.dat"), "w") as f:
        f.write(f"{len(at['t'])} {int(at.get('fourier_modes', 1))}\n")
        f.writelines(f"{t:.17g} {T:.17g}\n" for t, T in zip(at["t"], at["T"]))
    with open(os.path.join(case_dir, "svzerod_3Dcoupling.json"), "w") as f:
        json.dump(coupling_config(cfg, scale, n_substeps), f, indent=2)
    # GMRES: fsils' BiCGStab leaves out the coupled-BC resistance dP/dQ (the
    # face term of the tangent), which makes Newton a fixed-point iteration
    # that diverges for stiff circulations (closed valves). At most
    # max_iter restarts x krylov iterations. A tighter relative tolerance
    # over-solves near Newton convergence (round-off) and spoils the update.
    ls = dict(type="GMRES", algebra="fsils", preconditioner="fsils", tol=1e-6,
              abs_tol=1e-20, max_iter=4, krylov=300)
    ls.update(linear_solver or {})
    p0 = float(cfg.get("initial_condition", {}).get("pressure_all", 0.0))
    case = dict(source_0d=os.path.basename(model_0d), block=cyl["vessel_name"],
                Ri=v["Ri"], Re=v["Re"], length=v["length"], mesh_length=Lm, material=mat,
                fiber_files=[os.path.relpath(p, case_dir) for p in paths[:2]],
                fiber_location="points", robin_stiffness=0.0, formulation="struct-p2",
                order=fm["order"], sector=sector, nr=nr, ntheta=mesh["ntheta"], nz=mesh["nz"],
                n_nodes=len(mesh["points"]), n_cells=len(mesh["cells"]),
                faces=sorted(mesh["faces"]), period=period, cycles=cycles,
                n_steps=cycles * (n_pts - 1), dt=period / (n_pts - 1), density=density,
                newton_tol=newton_tol, linear_solver=ls, scale=scale,
                active_tension_file="active_tension.dat",
                coupling=dict(block=BLOCK, config="svzerod_3Dcoupling.json",
                              library=LIB_IN_IMAGE, initial_pressure=p0),
                # post-processing of the pressure ramp does not apply
                T_ramp=float("inf"), P_max=0.0,
                **({"end_damping": end_damping * mat["mu"] / v["Re"] * period / (n_pts - 1)}
                   if helical else {}))
    with open(os.path.join(case_dir, "solver.xml"), "w") as f:
        f.write(C3.solver_xml(case))
    with open(os.path.join(case_dir, "case.json"), "w") as f:
        json.dump(case, f, indent=2)
    return case


def run(case_dir, np_=1, max_time=None, echo_every=1):
    """Run the coupled case in Docker with this repository mounted at /svzerod."""
    lib = os.path.join(REPO, "build-docker", "src", "interface", "libsvzero_interface.so")
    if not os.path.exists(lib):
        raise FileNotFoundError(f"{lib} missing: build svzero_interface in Docker (see module doc)")
    return C3.run(case_dir, np_, mounts={REPO: "/svzerod"}, max_time=max_time,
                  echo_every=echo_every)


def read_svzerod_data(path, dt):
    """Time and named variables of svMultiPhysics' results/svZeroD_data.

    Each row holds the 0D state at the end of a time step but is labeled with
    the step's start time (svZeroD_subroutines.cpp: write_svZeroD_solution
    is called with svZeroDTime before it is advanced by dt), so dt is added.
    """
    with open(path) as f:
        header = f.readline().split()
        names = header[1:]
        rows = [list(map(float, ln.split())) for ln in f if ln.strip()]
    data = np.array(rows)
    return data[:, 0] + dt, {n: data[:, 1 + i] for i, n in enumerate(names)}


def post(case_dir):
    """Cavity volume, pressure, valve flows and wall kinematics per time step.

    Volume and flows are those of the full cylinder (3D volume times the
    model fraction, 0D flows times the impedance scale). With helical fibers
    (full cylinder, free ends) the wall kinematics and strains are those of
    the mid-length section (middle element layer, see cylinder_3d.post).
    Written to CASE_DIR/coupled.json.
    """
    import glob
    with open(os.path.join(case_dir, "case.json")) as f:
        case = json.load(f)
    Lm = case["mesh_length"]
    mesh = C3.mesh_cylinder(case["Ri"], case["Re"], Lm, case["nr"], case["sector"], case["order"],
                            case["ntheta"], case["nz"])
    X = mesh["points"]
    conn = mesh["faces"]["endo"][0]
    dR = (case["Re"] - case["Ri"]) / (case["nr"] * case["order"])
    levels = np.rint((np.hypot(X[:, 0], X[:, 1]) - case["Ri"]) / dR).astype(int)
    files = sorted(glob.glob(os.path.join(case_dir, "results", "result_*.vtu")),
                   key=lambda fn: int(os.path.basename(fn)[7:-4]))  # step number
    steps = [int(os.path.basename(fn)[7:-4]) for fn in files]
    out = dict(R=(case["Ri"] + dR * np.arange(levels.max() + 1)).tolist(),
               t=[0.0], V=[C3.lumen_volume(X, conn, case["sector"]) * case["length"] / Lm],
               rho=[[0.0] * (levels.max() + 1)], beta=[0.0])
    m = case["material"]
    fiber = C3.cylinder_fibers_at_nodes(X, case["Ri"], case["Re"], m.get("alpha_endo", 0.0),
                                        m.get("alpha_epi", 0.0))["fiber"]
    section, use = {}, None
    if case["sector"] == 360.0:  # mid-length section of the free-ended tube
        section = dict(z0=0.5 * Lm, dz=0.5 * Lm / case["nz"])
        use = np.abs(X[:, 2] - section["z0"]) <= section["dz"] * (1.0 + 1e-9)
        out["section"] = "mid-length"

    def fields(u):  # through-wall profiles of strain, Jacobian and fiber stretch
        kf = C3.kinematic_fields(X, mesh["cells"], u, case["order"], fiber)
        return {k: v.tolist() for k, v in
                C3.level_means(levels, kf, C3.KINEMATIC_PROFILES, use).items()}
    out["fields"] = {k: [v] for k, v in fields(np.zeros_like(X)).items()}
    for n, fn in zip(steps, files):
        pts, pdata = C3._read_vtu(fn)
        x = X + pdata["Displacement"][:, :3]
        k = C3.kinematics(X, x, levels, **section)
        out["t"].append(n * case["dt"])
        out["V"].append(C3.lumen_volume(x, conn, case["sector"]) * case["length"] / Lm)
        out["rho"].append(k["rho"].tolist())
        out["beta"].append(k["beta"])
        for key, v in fields(x - X).items():
            out["fields"][key].append(v)
    t0d, d = read_svzerod_data(os.path.join(case_dir, "results", "svZeroD_data"), case["dt"])
    s = case["scale"]
    out.update(t_0d=t0d.tolist(), P=d[f"pressure:{BLOCK}:cavity_3d"].tolist(),
               Q_in=(s * d["flow:inlet_valve:cavity_in"]).tolist(),
               Q_out=(s * d["flow:cavity_out:outlet_valve"]).tolist(),
               Q_3d=(s * d[f"flow:{BLOCK}:cavity_3d"]).tolist())
    with open(os.path.join(case_dir, "coupled.json"), "w") as f:
        json.dump(out, f)
    return out


def run_0d(model_0d, num_elements=None):
    """0D reference: the ChamberCylinder in the same circulation."""
    import svzerod
    with open(model_0d) as f:
        cfg = json.load(f)
    cyl = find_cylinder(cfg)
    if num_elements:
        cyl["zero_d_element_values"]["num_elements"] = num_elements
    name = cyl["vessel_name"]
    r = svzerod.simulate(cfg)

    def series(var):
        s = r[r["name"] == var].sort_values("time")
        return s["time"].to_numpy(), s["y"].to_numpy()
    t, V = series(f"volume:{name}")
    pname = next(n for n in r["name"].unique() if n.startswith("pressure:") and name in n)
    ne = int(cyl["zero_d_element_values"]["num_elements"])
    v = cyl["zero_d_element_values"]
    prof = {f: np.stack([series(f"{f}_{A}:{name}")[1] for A in range(ne + 1)], axis=1)
            for f in ("rho", "phi", "eta")}
    inlet = next(n for n in r["name"].unique() if n.startswith("flow:") and n.endswith(":" + name))
    outlet = next(n for n in r["name"].unique() if n.startswith(f"flow:{name}:"))
    return dict(t=t, V=V, P=series(pname)[1], Q_in=series(inlet)[1], Q_out=series(outlet)[1],
                beta=series(f"beta:{name}")[1], R=np.linspace(v["Ri"], v["Re"], ne + 1), **prof)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate")
    g.add_argument("model_0d")
    g.add_argument("case_dir")
    g.add_argument("--nr", type=int, required=True)
    g.add_argument("--cycles", type=int)
    r = sub.add_parser("run")
    r.add_argument("case_dir")
    r.add_argument("--np", type=int, default=1)
    r.add_argument("--max-time", type=float)
    a = ap.parse_args()
    if a.cmd == "generate":
        c = generate(a.model_0d, a.case_dir, a.nr, a.cycles)
        print(f"{a.case_dir}: {c['n_cells']} tets, {c['n_steps']} steps of {c['dt']:.4g} s, "
              f"flow scale {c['scale']:.0f}")
    elif a.cmd == "run":
        run(a.case_dir, a.np, a.max_time)


if __name__ == "__main__":
    main()
