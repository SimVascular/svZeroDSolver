"""Mesh convergence of the 0D block itself, and the quality of the reference.

Runs the 0D input of a study case at the radial resolutions --ne (cached in
CASE/output like the studies) and compares, at the compared time points, each
direct solution and each Richardson extrapolation (ne, 2 ne) with the best
estimate of the converged 0D solution: the Richardson extrapolation of the two
finest levels. Max relative differences of V, rho(R) and beta (passive, over
the pressure ramp) or V, P (relative to its range), Q_in (relative to its
peak), rho(R) and beta (active, --active: first cycle from 0.02 s). Profiles
are compared at the nodes of the coarser solution (the 0D meshes are nested).
Writes CASE/output/zero_d_convergence.json.

Usage: python zero_d_convergence.py CASE_DIR --ne 16 32 64 128 256 [--active]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cylinder_3d as C3  # noqa: E402
import run_study as RS  # noqa: E402


def richardson(c, f, keys):
    """Second-order Richardson extrapolation of the 0D solutions c (ne) and
    f (2 ne), with the profiles at the nodes of c."""
    r = dict(c)
    for k in keys:
        r[k] = (4.0 * f[k] - c[k]) / 3.0
    for k in ("rho", "phi", "eta"):
        if k in c:
            r[k] = (4.0 * f[k][:, ::2] - c[k]) / 3.0
    return r


def differences(sol, ref, m, active=False, twisting=True):
    """Max relative differences of sol vs ref at the time points m (beta
    only when the reference twists, RS.has_twist)."""
    R = sol["R"] if len(sol["R"]) <= len(ref["R"]) else ref["R"]
    on = lambda s: np.array([np.interp(R, s["R"], row) for row in s["rho"][m]])  # noqa: E731
    rho, rho_ref = on(sol), on(ref)
    e = dict(V=np.max(np.abs(sol["V"][m] / ref["V"][m] - 1.0)))
    if active:
        e["P"] = np.max(np.abs(sol["P"][m] - ref["P"][m])) / np.ptp(ref["P"][m])
        e["Q_in"] = np.max(np.abs(sol["Q_in"][m] - ref["Q_in"][m])) / np.max(np.abs(ref["Q_in"][m]))
    e["rho"] = np.max(np.abs(rho - rho_ref)) / np.max(np.abs(rho_ref))
    b = np.max(np.abs(ref["beta"][m]))
    if twisting and b > 0.0:
        e["beta"] = np.max(np.abs(sol["beta"][m] - ref["beta"][m])) / b
    return {k: float(x) for k, x in e.items()}


def study(root, sols, m, keys, L, active=False):
    """Differences of all levels and extrapolations to the finest extrapolation."""
    nes = sorted(sols)
    best = richardson(sols[nes[-2]], sols[nes[-1]], keys)
    tw = RS.has_twist(best, L)
    rows = {}
    print(f"{os.path.basename(root)}: 0D vs Richardson {nes[-2]}/{nes[-1]} "
          "(max relative difference)")
    for n in nes:
        rows[f"ne{n}"] = differences(sols[n], best, m, active, tw)
        print(f"  ne = {n:4d}: " + ", ".join(f"{k} {x:.2e}" for k, x in rows[f'ne{n}'].items()))
    for n in nes[:-2]:
        if 2 * n in sols:
            name = f"richardson{n}/{2 * n}"
            rows[name] = differences(richardson(sols[n], sols[2 * n], keys), best, m, active, tw)
            print(f"  Richardson {n}/{2 * n}: " + ", ".join(f"{k} {x:.2e}" for k, x in rows[name].items()))
    with open(os.path.join(root, "output", "zero_d_convergence.json"), "w") as f:
        json.dump(dict(case=os.path.basename(root), reference=f"richardson{nes[-2]}/{nes[-1]}",
                       differences=rows), f, indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("--ne", type=int, nargs="+", default=[16, 32, 64, 128, 256])
    ap.add_argument("--active", action="store_true", help="active case (run_active_study.py)")
    a = ap.parse_args()
    sys.stdout.reconfigure(line_buffering=True)
    root = os.path.abspath(a.case)
    model = os.path.join(root, "chamber_cylinder_0d.json")
    out_dir = os.path.join(root, "output")
    os.makedirs(out_dir, exist_ok=True)
    with open(model) as f:
        L = next(v for v in json.load(f)["vessels"]
                 if v.get("zero_d_element_type") == "ChamberCylinder")["zero_d_element_values"]["length"]
    if a.active:
        import run_active_study as RA
        sols = {}
        for n in sorted(a.ne):
            sols[n] = RA.zero_d(model, out_dir, n)  # cached in out_dir
            print(f"0D ne = {n} done", flush=True)
        t = sols[max(sols)]["t"]
        m = (t >= RA.STARTUP) & (t <= 1.0 + 1e-9)  # first cycle
        return study(root, sols, m, ("V", "P", "Q_in", "Q_out", "beta"), L, active=True)
    m0 = C3.read_0d(model)
    T, _ = C3.pressure_ramp(m0["t"], m0["P"])
    sp = m0["config"]["simulation_parameters"]
    t_all = np.linspace(0.0, sp["cardiac_period"], int(sp["number_of_time_pts_per_cardiac_cycle"]))
    t0 = t_all[t_all <= T * (1 + 1e-12)]
    sols = {ne: RS.restrict(RS.run_0d(m0["config"], m0["name"], ne, out_dir), t0) for ne in sorted(a.ne)}
    study(root, sols, np.ones(len(t0), bool), ("V", "beta"), L)


if __name__ == "__main__":
    main()
