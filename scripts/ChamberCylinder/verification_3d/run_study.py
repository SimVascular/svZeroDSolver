# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Passive inflation of a fixed-length cylinder: 0D vs. 3D mesh convergence.

CASE is a folder with the 0D input chamber_cylinder_0d.json, which defines the
problem. For each 3D mesh level (--nr, elements through the wall) the script
writes the svMultiPhysics case into CASE/runs/, runs it if --run is given
(Docker image simvascular/solver by default, or --exe), and compares over the
pressure ramp the kinematic quantities of the 0D block:

  V       cavity volume                        |V - V_ref| / V_ref
  rho     radial displacement rho(R)           max |rho - rho_ref| / max |rho_ref|
  h       wall-thickness change rho(Re)-rho(Ri) relative to its maximum
  beta    twist per unit length                relative to max |beta_ref|
  phi     transmural rotation phi(R)-phi(Ri)   relative to the twist max |beta_ref| L
  eta     axial displacement eta(R)-eta(Ri)    relative to max |rho_ref|

The profiles are compared at the radial node levels of the 3D mesh, which are
nodes of the 0D reference. The reference is the Richardson extrapolation of the
0D block at num_elements and 2 num_elements (second order), or the closed form
of the incompressible neo-Hookean tube when the wall is isotropic and
incompressible. Each error is the maximum over the 0D time points of the ramp;
the 3D model takes --substeps load steps per 0D step.

Run from the repository root after building the Python module:
  PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/run_study.py \\
      scripts/ChamberCylinder/verification_3d/passive_isotropic --nr 1 2 4 [--run] [--np 4]
Results: CASE/runs/<formulation>_nr<N>/ and CASE/output/ (summary.json, figure).
Pressures are reported in mmHg; the solver inputs are SI.
"""
import argparse
import copy
import glob
import hashlib
import json
import math
import os
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import cylinder_3d as C3  # noqa: E402
from ChamberCylinder_verification import fixed_length_exact  # noqa: E402

PROFILES = ("rho", "phi", "eta")
QUANTITIES = [  # key, label, plot marker
    ("V", "volume", "s-"), ("rho", "radial displacement", "o-"),
    ("h", "wall-thickness change", "D-"), ("beta", "twist", "^-"),
    ("phi", "transmural rotation", "v-"), ("eta", "axial displacement", "x-")]


def run_0d(config, name, num_elements=None, cache_dir=None):
    """0D pressure, volume, twist and nodal profiles rho, phi, eta over time.

    Results are cached in cache_dir under a hash of the 0D input.
    """
    cfg = copy.deepcopy(config)
    blk = next(v for v in cfg["vessels"] if v["vessel_name"] == name)
    vals = blk["zero_d_element_values"]
    if num_elements:
        vals["num_elements"] = num_elements
    ne = int(vals.get("num_elements", 10))
    cache = None
    if cache_dir:
        key = hashlib.sha1(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16]
        cache = os.path.join(cache_dir, f"zero_d_{ne}_{key}.npz")
        if os.path.exists(cache):
            with np.load(cache) as d:
                return {k: d[k] for k in d.files}
    import svzerod
    tic = time.time()
    out = svzerod.simulate(cfg)

    def series(var):
        d = out[out["name"] == var].sort_values("time")
        return d["time"].to_numpy(), d["y"].to_numpy()
    t, V = series(f"volume:{name}")
    pname = next(n for n in out["name"].unique() if n.startswith("pressure:") and name in n)
    sol = dict(t=t, P=series(pname)[1], V=V, beta=series(f"beta:{name}")[1],
               R=np.linspace(vals["Ri"], vals["Re"], ne + 1))
    for f in PROFILES:
        sol[f] = np.stack([series(f"{f}_{A}:{name}")[1] for A in range(ne + 1)], axis=1)
    if cache:
        os.makedirs(cache_dir, exist_ok=True)
        np.savez(cache, **sol)
    print(f"0D ne = {ne}: {time.time() - tic:.0f} s", flush=True)
    return sol


def exact_isotropic(t, P, Ri, Re, L, mu, n=257):
    """Closed form of the fixed-length incompressible neo-Hookean tube."""
    R = np.linspace(Ri, Re, n)
    ri = np.array([fixed_length_exact(p, Ri, Re, mu) for p in P])
    rho = np.sqrt(R[None, :] ** 2 + ri[:, None] ** 2 - Ri ** 2) - R[None, :]
    z = np.zeros_like(rho)
    return dict(t=t, P=P, V=math.pi * ri ** 2 * L, beta=np.zeros_like(P), R=R,
                rho=rho, phi=z, eta=z)


def richardson(coarse, fine):
    """Second-order Richardson extrapolation of two 0D solutions (ne, 2 ne) on
    the coarse nodes."""
    ext = dict(coarse)
    for k in ("V", "beta"):
        ext[k] = (4.0 * fine[k] - coarse[k]) / 3.0
    for k in PROFILES:
        ext[k] = (4.0 * fine[k][:, ::2] - coarse[k]) / 3.0
    return ext


def restrict(sol, t):
    """Solution at the times t (a subset of its time points)."""
    idx = [int(np.argmin(np.abs(sol["t"] - ti))) for ti in t]
    if not np.allclose(sol["t"][idx], t, rtol=0, atol=1e-9 * max(t)):
        raise ValueError("the compared solutions have different time points")
    pick = lambda v: {k: pick(w) for k, w in v.items()} if isinstance(v, dict) \
        else np.asarray(v)[idx]  # noqa: E731
    return {k: (v if k == "R" else pick(v)) for k, v in sol.items()}


def tolist(x):
    """JSON-ready copy of an array or a dict of arrays."""
    return {k: tolist(v) for k, v in x.items()} if isinstance(x, dict) else np.asarray(x).tolist()


def has_twist(ref, L):
    """True when the reference twists beyond round-off (helical fibers)."""
    twist = np.max(np.abs(ref["beta"])) * L * np.max(ref["R"])  # end-face rotation
    return bool(twist > 1e-8 * np.max(np.abs(ref["rho"])))


def errors(sol, ref, L):
    """Max relative differences over the ramp of sol vs ref (same times)."""
    at = lambda f: np.array([np.interp(sol["R"], ref["R"], row) for row in ref[f]])  # noqa: E731
    rho, rho_ref = np.asarray(sol["rho"]), at("rho")
    dh, dh_ref = rho[:, -1] - rho[:, 0], rho_ref[:, -1] - rho_ref[:, 0]
    scale_rho = np.max(np.abs(rho_ref))
    e = dict(V=np.max(np.abs(np.asarray(sol["V"]) / ref["V"] - 1.0)),
             rho=np.max(np.abs(rho - rho_ref)) / scale_rho,
             h=np.max(np.abs(dh - dh_ref)) / np.max(np.abs(dh_ref)),
             eta=np.max(np.abs(np.asarray(sol["eta"]) - at("eta"))) / scale_rho)
    b_ref = np.max(np.abs(ref["beta"]))
    if has_twist(ref, L):  # helical fibers
        e["beta"] = np.max(np.abs(np.asarray(sol["beta"]) - ref["beta"])) / b_ref
        e["phi"] = np.max(np.abs(np.asarray(sol["phi"]) - at("phi"))) / (b_ref * L)
    return {k: float(v) for k, v in e.items()}


def table_row(label, e):
    return f"{label:<12s}" + "".join(f"{e[k]:>12.2e}" if k in e else f"{'-':>12s}"
                                     for k, _, _ in QUANTITIES)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case", help="folder with chamber_cylinder_0d.json")
    ap.add_argument("--nr", type=int, nargs="+", default=[1, 2, 4])
    ap.add_argument("--formulation", default="struct-p2", choices=sorted(C3.FORMULATIONS))
    ap.add_argument("--sector", type=float, help="default: 90, or 360 with helical fibers")
    ap.add_argument("--full-length", action="store_true",
                    help="mesh the full length (default: half, mid-length symmetry plane)")
    ap.add_argument("--run", action="store_true",
                    help="run svMultiPhysics where complete results of the same case are missing")
    ap.add_argument("--rerun", action="store_true", help="with --run: also rerun complete cases")
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--exe", help="local svmultiphysics executable instead of Docker")
    ap.add_argument("--max-time", type=float, default=1800.0,
                    help="wall-time limit per 3D run [s]")
    ap.add_argument("--substeps", type=int, default=1, help="3D load steps per 0D time step")
    ap.add_argument("--ls-type", default="BICG")
    ap.add_argument("--algebra", default="fsils")
    ap.add_argument("--preconditioner", default="fsils")
    a = ap.parse_args()
    sys.stdout.reconfigure(line_buffering=True)

    case_root = os.path.abspath(a.case)
    model_0d = os.path.join(case_root, "chamber_cylinder_0d.json")
    m0 = C3.read_0d(model_0d)
    v = m0["values"]
    mat = C3.material(v)
    fibers = mat["model"] in ("HGO", "HolzapfelOgden")
    T, Pmax = C3.pressure_ramp(m0["t"], m0["P"])
    Ri, Re, L = v["Ri"], v["Re"], v["length"]
    mmhg = C3.MMHG

    # time points of the 0D input on the ramp; the 3D case uses the same steps
    sp = m0["config"]["simulation_parameters"]
    t_all = np.linspace(0.0, sp["cardiac_period"], int(sp["number_of_time_pts_per_cardiac_cycle"]))
    t0 = t_all[t_all <= T * (1 + 1e-12)]
    out_dir = os.path.join(case_root, "output")
    zero_d = lambda ne=None: restrict(run_0d(m0["config"], m0["name"], ne, out_dir), t0)  # noqa: E731

    # 3D: generate and run each mesh level (stop at the first failure)
    ls = dict(type=a.ls_type, algebra=a.algebra, preconditioner=a.preconditioner)
    cases = {}
    for nr in a.nr:
        case_dir = os.path.join(case_root, "runs", f"{a.formulation}_nr{nr}")
        os.makedirs(case_dir, exist_ok=True)
        old = None
        if os.path.exists(os.path.join(case_dir, "case.json")):
            with open(os.path.join(case_dir, "case.json")) as f:
                old = json.load(f)
        case = C3.generate(model_0d, case_dir, nr, a.formulation, a.sector,
                           n_steps=a.substeps * (len(t0) - 1), linear_solver=ls,
                           half_length=not a.full_length)
        cases[nr] = (case_dir, case)
        done = len(glob.glob(os.path.join(case_dir, "results", "result_*.vtu"))) == case["n_steps"]
        if a.run and not a.rerun and done and old == json.loads(json.dumps(case)):
            print(f"3D nr = {nr}: complete results of the same case exist, not rerun")
        elif a.run:
            print(f"3D nr = {nr}: running ({case['n_cells']} tets, {case['n_steps']} steps)")
            shutil.rmtree(os.path.join(case_dir, "results"), ignore_errors=True)  # no stale steps
            try:
                wall = C3.run(case_dir, a.np, a.exe, max_time=a.max_time)
            except C3.SolverFailure as err:
                print(f"3D nr = {nr}: FAILED - {err}")
                break  # finer meshes would fail too
            with open(os.path.join(case_dir, "wall_time.txt"), "w") as f:
                f.write(f"{wall:.1f}\n")
            print(f"3D nr = {nr}: done in {wall:.0f} s")

    # 0D reference
    if mat["kappa"] > 0.0 or fibers:
        ne = int(v["num_elements"])
        ref, ref_name = richardson(zero_d(ne), zero_d(2 * ne)), f"0D Richardson ne = {ne}/{2 * ne}"
    else:
        ref, ref_name = exact_isotropic(t0, Pmax * t0 / T, Ri, Re, L, mat["mu"]), "exact"
    summary = dict(case=os.path.basename(case_root), material=mat, P_max_mmHg=Pmax / mmhg,
                   reference=dict(name=ref_name, **{k: np.asarray(x).tolist()
                                                    for k, x in ref.items()}),
                   zero_d={}, three_d={})
    twist = np.degrees(ref["beta"][-1] * L)
    print(f"0D reference ({ref_name}) at {Pmax / mmhg:.0f} mmHg: V = {ref['V'][-1] * 1e6:.4f} mL, "
          f"rho endo/epi = {ref['rho'][-1][0] * 1e3:.4f}/{ref['rho'][-1][-1] * 1e3:.4f} mm, "
          f"twist over the length = {twist:.4f} deg")

    # compare
    rows = []
    for nr, (case_dir, case) in cases.items():
        try:
            out = C3.post(case_dir)
        except FileNotFoundError:
            print(f"3D nr = {nr}: generated {case_dir} ({case['n_cells']} tets), no results yet")
            continue
        sol = {k: (x if isinstance(x, dict) else np.asarray(x)) for k, x in out.items()}
        if len(sol["t"]) != case["n_steps"] + 1:
            print(f"3D nr = {nr}: incomplete results ({len(sol['t']) - 1} of {case['n_steps']} steps)")
            continue
        sol = restrict(sol, t0)
        e = errors(sol, ref, L)
        rec = dict(h=(Re - Ri) / nr, n_cells=case["n_cells"], n_nodes=case["n_nodes"], errors=e,
                   non_homogeneity=float(np.max(sol["rho_spread"]) / np.max(np.abs(ref["rho"]))),
                   **{k: tolist(x) for k, x in sol.items()})
        wall = os.path.join(case_dir, "wall_time.txt")
        rec["wall_time"] = float(open(wall).read()) if os.path.exists(wall) else None
        summary["three_d"][f"nr{nr}"] = rec
        rows.append(table_row(f"3D nr = {nr}", e)
                    + f"   (non-homogeneity {rec['non_homogeneity']:.1e}"
                    + (f", {rec['wall_time']:.0f} s)" if rec["wall_time"] else ")"))
    for ne in (2, 4, 8, 16, 32):  # 0D mesh convergence against the same reference
        e = errors(zero_d(ne), ref, L)
        summary["zero_d"][f"ne{ne}"] = e
        rows.append(table_row(f"0D ne = {ne}", e))
    header = f"{'':<12s}" + "".join(f"{k:>12s}" for k, _, _ in QUANTITIES)
    print("max relative error over the ramp\n" + header + "\n" + "\n".join(rows))

    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    twisting = has_twist(ref, L)
    plot(summary, ref, ref_name, mat, L, twisting, a.formulation, out_dir)


def plot(summary, ref, ref_name, mat, L, twisting, formulation, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mmhg = C3.MMHG
    Pm = summary["P_max_mmHg"]
    levels = summary["three_d"]
    P0 = np.asarray(ref["P"]) / mmhg
    n = 4 if twisting else 3
    fig, axs = plt.subplots(1, n, figsize=(5.2 * n, 4.4))
    ax_v, ax_r, ax_c = axs[0], axs[1], axs[-1]

    def bands(a):
        a.axhspan(5, 12, color="0.9", lw=0)
        a.axhspan(100, 120, color="0.9", lw=0)
        a.text(0.02, 8.5 / Pm, "end-diastolic", transform=a.transAxes, va="center")
        a.text(0.02, 110 / Pm, "peak systolic", transform=a.transAxes, va="center")
        a.set(ylabel="cavity pressure [mmHg]", ylim=(0, Pm))

    # pressure-volume
    bands(ax_v)
    if mat["model"] == "neoHookean" and ref_name != "exact":  # compressible isotropic: show the limit
        Pe = np.linspace(0.0, Pm, 101)
        Ve = [math.pi * fixed_length_exact(p * mmhg, ref["R"][0], ref["R"][-1], mat["mu"]) ** 2 * L
              for p in Pe]
        ax_v.plot(np.array(Ve) * 1e6, Pe, "k:", lw=1.5, label="incompressible, exact")
    ax_v.plot(np.asarray(ref["V"]) * 1e6, P0, "k-", lw=2.5, label=ref_name)
    for i, (k, r) in enumerate(levels.items()):
        ax_v.plot(np.array(r["V"]) * 1e6, np.array(r["P"]) / mmhg, "--", c=f"C{i}", label=f"3D {k}")
    ax_v.set(xlabel="cavity volume [mL]", title="pressure-volume")
    ax_v.legend(loc="lower right")

    # radial displacement through the wall at end diastole and peak pressure
    R = np.asarray(ref["R"]) * 1e3
    for P_show in (10.0, Pm):
        j = int(np.argmin(np.abs(P0 - P_show)))
        ax_r.plot(R, np.asarray(ref["rho"])[j] * 1e3, "k-", lw=2.0)
        ax_r.annotate(f"{P0[j]:.0f} mmHg", (R[-1], ref["rho"][j][-1] * 1e3), ha="right",
                      va="bottom")
        for i, (k, r) in enumerate(levels.items()):
            jr = int(np.argmin(np.abs(np.array(r["P"]) / mmhg - P_show)))
            ax_r.plot(np.array(r["R"]) * 1e3, np.array(r["rho"][jr]) * 1e3, "o", ms=4, mfc="none",
                      c=f"C{i}", label=f"3D {k}" if P_show == Pm else None)
    ax_r.plot([], [], "k-", label=ref_name)
    ax_r.set(xlabel="reference radius R [mm]", ylabel="radial displacement rho [mm]",
             title="through-wall displacement")
    ax_r.legend()

    if twisting:  # twist
        ax_t = axs[2]
        bands(ax_t)
        ax_t.plot(np.degrees(np.asarray(ref["beta"]) * L), P0, "k-", lw=2.5, label=ref_name)
        for i, (k, r) in enumerate(levels.items()):
            ax_t.plot(np.degrees(np.array(r["beta"]) * L), np.array(r["P"]) / mmhg, "--",
                      c=f"C{i}", label=f"3D {k}")
        ax_t.set(xlabel="twist over the length [deg]", title="twist")
        ax_t.legend(loc="lower right")

    if levels:
        h = np.array([r["h"] for r in levels.values()]) * 1e3
        for key, lab, mk in QUANTITIES:
            err = [r["errors"].get(key) for r in levels.values()]
            if None not in err and max(err) > 0.0:
                ax_c.loglog(h, err, mk, label=lab)
        e0 = levels[next(iter(levels))]["errors"]["V"]
        ax_c.loglog(h, e0 * (h / h[0]) ** 2, "k:", label="slope 2")
        ax_c.loglog(h, e0 * (h / h[0]) ** 4, "k-.", label="slope 4")
    ax_c.set(xlabel="element size h = (Re - Ri)/nr [mm]", ylabel="max relative error vs 0D",
             title=f"3D mesh convergence ({formulation})")
    ax_c.grid(True, which="both", alpha=0.3)
    ax_c.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"convergence_{formulation}.png"), dpi=130)


if __name__ == "__main__":
    main()
