# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Active contraction in a circulation: 0D cylinder vs. coupled 3D cylinder.

CASE is a folder with the 0D input chamber_cylinder_0d.json: a ChamberCylinder
(prescribed active tension, see the block docs) in a circulation. For each 3D
mesh level (--nr) the script writes the coupled svMultiPhysics case
(coupled_3d.py) to CASE/runs/, runs it if --run is given (stopping at the
first failure; complete runs of the same case are reused), and compares over
the last cardiac cycle with the 0D block in the same circulation:

  V      cavity volume                   max |V - V_ref| / V_ref
  P      cavity pressure                 max |P - P_ref| / (max P_ref - min P_ref)
  rho    radial displacement rho(R)      max |rho - rho_ref| / max |rho_ref|
  SV     stroke volume                   |SV - SV_ref| / SV_ref

The 0D reference is the Richardson extrapolation of the block at 2 and 4
times num_elements of the input.

Run from the repository root after building the Python module and the Linux
svZeroDSolver interface library (see coupled_3d.py):
  PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/run_active_study.py \\
      scripts/ChamberCylinder/verification_3d/active_prescribed --nr 1 2 [--run] [--np 2]
"""
import argparse
import glob
import hashlib
import json
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import coupled_3d as CP  # noqa: E402

MMHG = 133.322387415


def cache_key(model_0d):
    """Hash of the 0D input without its description."""
    with open(model_0d) as f:
        cfg = json.load(f)
    cfg.pop("description", None)
    return hashlib.sha1(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16]


def reference(model_0d, cache_dir):
    """Richardson-extrapolated 0D solution (2 and 4 times num_elements), cached."""
    with open(model_0d) as f:
        ne = int(CP.find_cylinder(json.load(f))["zero_d_element_values"]["num_elements"])
    key = cache_key(model_0d)
    sols = []
    for n in (2 * ne, 4 * ne):
        path = os.path.join(cache_dir, f"zero_d_ne{n}_{key}.npz")
        if os.path.exists(path):
            with np.load(path) as d:
                sols.append({k: d[k] for k in d.files})
            continue
        s = CP.run_0d(model_0d, n)
        os.makedirs(cache_dir, exist_ok=True)
        np.savez(path, **s)
        sols.append(s)
    c, f = sols
    ref = dict(c)
    for k in ("V", "P", "Q_in", "Q_out"):
        ref[k] = (4.0 * f[k] - c[k]) / 3.0
    ref["rho"] = (4.0 * f["rho"][:, ::2] - c["rho"]) / 3.0
    return ref, f"0D Richardson ne = {2 * ne}/{4 * ne}"


def last_cycle(t, period, cycles):
    return t >= (cycles - 1) * period - 1e-9


def errors(sol, ref, period, cycles):
    """Errors over the last cycle (3D and 0D on the same time points)."""
    t3 = np.asarray(sol["t"])
    idx = [int(np.argmin(np.abs(ref["t"] - t))) for t in t3]
    if not np.allclose(ref["t"][idx], t3, rtol=0, atol=1e-9):
        raise ValueError("0D and 3D time points differ")
    m = last_cycle(t3, period, cycles)
    V, Vr = np.asarray(sol["V"])[m], ref["V"][idx][m]
    P = np.interp(t3, sol["t_0d"], sol["P"])[m]
    Pr = ref["P"][idx][m]
    rho = np.asarray(sol["rho"])[m]
    rr = np.array([np.interp(sol["R"], ref["R"], row) for row in ref["rho"][idx][m]])
    return dict(V=float(np.max(np.abs(V / Vr - 1.0))),
                P=float(np.max(np.abs(P - Pr)) / (Pr.max() - Pr.min())),
                rho=float(np.max(np.abs(rho - rr)) / np.max(np.abs(rr))),
                SV=float(abs((V.max() - V.min()) / (Vr.max() - Vr.min()) - 1.0)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case")
    ap.add_argument("--nr", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--rerun", action="store_true")
    ap.add_argument("--np", type=int, default=2)
    ap.add_argument("--max-time", type=float, default=7200.0)
    a = ap.parse_args()
    sys.stdout.reconfigure(line_buffering=True)
    root = os.path.abspath(a.case)
    model = os.path.join(root, "chamber_cylinder_0d.json")
    with open(model) as f:
        sp = json.load(f)["simulation_parameters"]
    period, cycles = float(sp["cardiac_period"]), int(sp["number_of_cardiac_cycles"])
    out_dir = os.path.join(root, "output")
    sols = {}
    for nr in a.nr:
        d = os.path.join(root, "runs", f"struct-p2_nr{nr}")
        old = None
        if os.path.exists(os.path.join(d, "case.json")):
            with open(os.path.join(d, "case.json")) as f:
                old = json.load(f)
        case = CP.generate(model, d, nr)
        done = len(glob.glob(os.path.join(d, "results", "result_*.vtu"))) == case["n_steps"]
        if a.run and (a.rerun or not done or old != json.loads(json.dumps(case))):
            print(f"3D nr = {nr}: running ({case['n_cells']} tets, {case['n_steps']} steps)")
            shutil.rmtree(os.path.join(d, "results"), ignore_errors=True)  # no stale steps
            try:
                wall = CP.run(d, a.np, a.max_time, echo_every=50)
            except CP.C3.SolverFailure as err:
                print(f"3D nr = {nr}: FAILED - {err}")
                break
            with open(os.path.join(d, "wall_time.txt"), "w") as f:
                f.write(f"{wall:.1f}\n")
            print(f"3D nr = {nr}: done in {wall:.0f} s")
        if len(glob.glob(os.path.join(d, "results", "result_*.vtu"))) != case["n_steps"]:
            print(f"3D nr = {nr}: no complete results")
            continue
        sols[nr] = CP.post(d)
    ref, ref_name = reference(model, out_dir)
    m = last_cycle(ref["t"], period, cycles)
    print(f"{ref_name}, last cycle: V {ref['V'][m].min() * 1e6:.2f}-{ref['V'][m].max() * 1e6:.2f} mL,"
          f" P {ref['P'][m].min() / MMHG:.1f}-{ref['P'][m].max() / MMHG:.1f} mmHg")
    summary = dict(reference=ref_name, three_d={})
    print(f"{'':10s}{'V':>11s}{'P':>11s}{'rho':>11s}{'SV':>11s}")
    for nr, sol in sols.items():
        e = errors(sol, ref, period, cycles)
        summary["three_d"][f"nr{nr}"] = dict(errors=e)
        print(f"3D nr = {nr:<3d}" + "".join(f"{e[k]:11.2e}" for k in ("V", "P", "rho", "SV")))
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    plot(sols, ref, ref_name, period, cycles, out_dir)


def plot(sols, ref, ref_name, period, cycles, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15.6, 4.4))
    m = last_cycle(ref["t"], period, cycles)
    t0 = (cycles - 1) * period
    ax[0].plot(ref["V"][m] * 1e6, ref["P"][m] / MMHG, "k-", lw=2.5, label=ref_name)
    ax[1].plot(ref["t"][m] - t0, ref["V"][m] * 1e6, "k-", lw=2.5, label=ref_name)
    ax[2].plot(ref["t"][m] - t0, ref["P"][m] / MMHG, "k-", lw=2.5, label=ref_name)
    for i, (nr, s) in enumerate(sols.items()):
        t = np.asarray(s["t"])
        k = last_cycle(t, period, cycles)
        P = np.interp(t, s["t_0d"], s["P"])
        lab = f"3D nr{nr}"
        ax[0].plot(np.asarray(s["V"])[k] * 1e6, P[k] / MMHG, "--", c=f"C{i}", label=lab)
        ax[1].plot(t[k] - t0, np.asarray(s["V"])[k] * 1e6, "--", c=f"C{i}", label=lab)
        ax[2].plot(t[k] - t0, P[k] / MMHG, "--", c=f"C{i}", label=lab)
    ax[0].set(xlabel="cavity volume [mL]", ylabel="cavity pressure [mmHg]", title="PV loop")
    ax[1].set(xlabel="time in cycle [s]", ylabel="cavity volume [mL]", title="volume")
    ax[2].set(xlabel="time in cycle [s]", ylabel="cavity pressure [mmHg]", title="pressure")
    for x in ax:
        x.legend(fontsize=8)
        x.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "pv_loop.png"), dpi=130)


if __name__ == "__main__":
    main()
