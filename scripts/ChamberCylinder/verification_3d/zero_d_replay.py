"""Diagnostic: the 0D block loaded with the cavity pressure of a coupled 3D run.

In the coupled comparison (run_active_study.py) the circulation feeds every
difference of the 3D tube back into the pressure: with helical fibers, the end
layers of the free-ended 3D tube change the whole-tube volume, hence the
pressure, hence the wall state everywhere. Here the 0D block alone (the same
wall and active tension, closed outlet) is loaded with the 3D cavity pressure
P3(t) of results/svZeroD_data, so both models are compared at the same
pressure: the whole-tube volume difference is the effect of the free ends,
and the mid-length differences are what remains of the wall mechanics. The 3D
model is not changed; this is post-processing of its results.

Writes CASE/output/replay.json (errors over the first cycle, from 0.02 s on).

Usage: python zero_d_replay.py CASE_DIR --nr 1 2 [--ne N]
"""
import argparse
import copy
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import coupled_3d as CP  # noqa: E402
import run_active_study as RA  # noqa: E402


def replay_config(cfg, t, P):
    """0D input: the cylinder of cfg between a pressure P(t) and a closed outlet."""
    c = copy.deepcopy(cfg)
    cyl = CP.find_cylinder(c)
    cyl["boundary_conditions"] = {"inlet": "CAVITY", "outlet": "CLOSED"}
    c["vessels"], c["valves"], c["junctions"] = [cyl], [], []
    dt = t[1] - t[0]
    # one sample past the end: the time-dependent BC is periodic over its samples
    tb, Pb = np.append(t, t[-1] + dt), np.append(P, P[-1])
    c["boundary_conditions"] = [
        {"bc_name": "CAVITY", "bc_type": "PRESSURE", "bc_values": {"t": tb.tolist(), "P": Pb.tolist()}},
        {"bc_name": "CLOSED", "bc_type": "FLOW", "bc_values": {"t": [0.0, float(tb[-1])], "Q": [0.0, 0.0]}}]
    sp = c["simulation_parameters"]
    sp.update(number_of_cardiac_cycles=1, cardiac_period=float(tb[-1]),
              number_of_time_pts_per_cardiac_cycle=len(tb))
    ic = c.setdefault("initial_condition", {})
    ic["pressure_all"] = float(P[0])
    return c


def replay_0d(model, ne=None):
    """The replay solution: the direct 0D solution with ne radial elements, or
    (default) the Richardson extrapolation of 2 and 4 times num_elements, as
    the reference of run_active_study.py."""
    if ne:
        return CP.run_0d(model, ne)
    with open(model) as f:
        n = int(CP.find_cylinder(json.load(f))["zero_d_element_values"]["num_elements"])
    c, f = CP.run_0d(model, 2 * n), CP.run_0d(model, 4 * n)
    z = dict(c)
    z["V"], z["beta"] = ((4.0 * f[k] - c[k]) / 3.0 for k in ("V", "beta"))
    z["rho"] = (4.0 * f["rho"][:, ::2] - c["rho"]) / 3.0
    return z


def compare(sol, z, L):
    """Max relative differences over the first cycle (from RA.STARTUP)."""
    t3 = np.asarray(sol["t"])
    m = (t3 >= RA.STARTUP) & (t3 <= 1.0 + 1e-9)
    idx = [int(np.argmin(np.abs(z["t"] - x))) for x in t3[m]]
    if not np.allclose(z["t"][idx], t3[m], rtol=0, atol=1e-9):
        raise ValueError("0D and 3D time points differ")
    V3, V0 = np.asarray(sol["V"])[m], z["V"][idx]
    rho3 = np.asarray(sol["rho"])[m]
    rho0 = np.array([np.interp(sol["R"], z["R"], row) for row in z["rho"][idx]])
    Ri = sol["R"][0]
    Vmid = np.pi * (Ri + rho3[:, 0]) ** 2 * L  # the 0D volume formula at mid-length
    sv = lambda V: V.max() - V.min()  # noqa: E731
    e = dict(V=np.max(np.abs(V3 / V0 - 1.0)), SV=abs(sv(V3) / sv(V0) - 1.0),
             V_mid=np.max(np.abs(Vmid / V0 - 1.0)), SV_mid=abs(sv(Vmid) / sv(V0) - 1.0),
             rho=np.max(np.abs(rho3 - rho0)) / np.max(np.abs(rho0)))
    b0 = z["beta"][idx]
    if np.max(np.abs(b0)) > 1e-12:
        e["beta"] = np.max(np.abs(np.asarray(sol["beta"])[m] - b0)) / np.max(np.abs(b0))
    return {k: float(x) for k, x in e.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("--nr", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--ne", type=int, help="direct 0D solution with this many radial elements "
                                           "(default: Richardson extrapolation, as the reference)")
    a = ap.parse_args()
    root = os.path.abspath(a.case)
    with open(os.path.join(root, "chamber_cylinder_0d.json")) as f:
        cfg = json.load(f)
    L = CP.find_cylinder(cfg)["zero_d_element_values"]["length"]
    out_dir = os.path.join(root, "output")
    os.makedirs(out_dir, exist_ok=True)
    n = int(CP.find_cylinder(cfg)["zero_d_element_values"]["num_elements"])
    result = dict(case=os.path.basename(root), levels={},
                  reference=f"0D ne = {a.ne}" if a.ne else f"0D Richardson ne = {2 * n}/{4 * n}")
    for nr in a.nr:
        d = os.path.join(root, "runs", f"struct-p2_nr{nr}")
        path = os.path.join(d, "coupled.json")
        if not os.path.exists(path):
            print(f"nr = {nr}: no coupled.json (run run_active_study.py first)")
            continue
        with open(path) as f:
            sol = json.load(f)
        P0 = json.load(open(os.path.join(d, "case.json")))["coupling"]["initial_pressure"]
        t = np.concatenate([[0.0], sol["t_0d"]])
        P = np.concatenate([[P0], sol["P"]])
        model = os.path.join(out_dir, f"replay_0d_nr{nr}.json")
        with open(model, "w") as f:
            json.dump(replay_config(cfg, t, P), f, indent=1)
        z = replay_0d(model, a.ne)
        result["levels"][f"nr{nr}"] = e = compare(sol, z, L)
        print(f"nr = {nr}: " + ", ".join(f"{k} {x:.2e}" for k, x in e.items()))
    with open(os.path.join(out_dir, "replay.json"), "w") as f:
        json.dump(result, f, indent=1)


if __name__ == "__main__":
    main()
