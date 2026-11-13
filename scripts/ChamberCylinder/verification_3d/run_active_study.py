# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Active contraction in a circulation: 0D cylinder vs. coupled 3D cylinder.

CASE is a folder with the 0D input chamber_cylinder_0d.json: a ChamberCylinder
(prescribed active tension, see the block docs) in a circulation. For each 3D
mesh level (--nr) the script writes the coupled svMultiPhysics case
(coupled_3d.py) to CASE/runs/, runs it if --run is given (stopping at the
first failure; complete runs of the same case are reused), and compares over
the last cardiac cycle with the 0D block in the same circulation:

  V       cavity volume                  max |V - V_ref| / V_ref
  P       cavity pressure                max |P - P_ref| / (max P_ref - min P_ref)
  rho     radial displacement rho(R)     max |rho - rho_ref| / max |rho_ref|
  SV      stroke volume                  |SV - SV_ref| / SV_ref
  Q_in    inlet valve flow               max |Q - Q_ref| / max |Q_ref|
  Q_out   outlet valve flow              max |Q - Q_ref| / max |Q_ref|
  W       stroke work, PV-loop area      |W - W_ref| / W_ref
  EF      ejection fraction SV / V_max   |EF - EF_ref| / EF_ref
  dPdt+   peak dP/dt (contraction)       |dPdt - dPdt_ref| / dPdt_ref
  dPdt-   minimum dP/dt (relaxation)     |dPdt - dPdt_ref| / dPdt_ref
  valves  opening and closing times      max |t - t_ref| [ms]

The valve times are the crossings of 2% of the peak flow of each valve in the
reference. Through the wall, the Green-Lagrange strains E_RR, E_TT, E_ZZ, the
Jacobian J and the fiber stretch lambda_f are compared as in run_study.py
(level means of the 3D nodal fields vs the 0D profiles, J as element means).

The 0D reference is the Richardson extrapolation of the block at 2 and 4
times num_elements of the input.

Run from the repository root after building the Python module and the Linux
svZeroDSolver interface library (see coupled_3d.py):
  PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/run_active_study.py \\
      scripts/ChamberCylinder/verification_3d/active_prescribed --nr 1 2 [--run] [--np 2]
Results: CASE/output/ (summary.json, comparison_cycle<k>.pdf).
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
import run_study as RS  # noqa: E402

MMHG = 133.322387415


def cache_key(model_0d):
    """Hash of the 0D input without its description."""
    with open(model_0d) as f:
        cfg = json.load(f)
    cfg.pop("description", None)
    return hashlib.sha1(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16]


def zero_d(model_0d, cache_dir, n):
    """0D solution with n radial elements, cached."""
    path = os.path.join(cache_dir, f"zero_d_ne{n}_{cache_key(model_0d)}.npz")
    if os.path.exists(path):
        with np.load(path) as d:
            s = {k: d[k] for k in d.files}
        if all(k in s for k in ("beta", "phi", "eta")):  # complete cache
            return s
    s = CP.run_0d(model_0d, n)
    os.makedirs(cache_dir, exist_ok=True)
    np.savez(path, **s)
    return s


def reference(model_0d, cache_dir, direct_ne=None):
    """Richardson-extrapolated 0D solution (2 and 4 times num_elements), or
    the direct solution with direct_ne elements."""
    if direct_ne:
        return zero_d(model_0d, cache_dir, direct_ne), f"0D ne = {direct_ne}"
    with open(model_0d) as f:
        ne = int(CP.find_cylinder(json.load(f))["zero_d_element_values"]["num_elements"])
    c, f = (zero_d(model_0d, cache_dir, n) for n in (2 * ne, 4 * ne))
    ref = dict(c)
    for k in ("V", "P", "Q_in", "Q_out", "beta"):
        ref[k] = (4.0 * f[k] - c[k]) / 3.0
    for k in ("rho", "phi", "eta"):
        ref[k] = (4.0 * f[k][:, ::2] - c[k]) / 3.0
    return ref, f"0D Richardson ne = {2 * ne}/{4 * ne}"


STARTUP = 0.02  # s, excluded from the comparison of the first cycle


def cycle_mask(t, period, k):
    """Time points of cycle k (1-based), both ends included."""
    return (t >= (k - 1) * period - 1e-9) & (t <= k * period + 1e-9)


# Linear solvers of the coupled cases (--ls). fsils' BiCGStab drops the
# coupled-BC tangent, so both use GMRES; fsils' diagonal preconditioner does
# not scale beyond nr = 2, Trilinos' algebraic multigrid does (nr = 4: ~180
# iterations per solve).
LINEAR_SOLVERS = {
    "fsils": None,  # coupled_3d default: GMRES(300) x 4, tol 1e-6
    "trilinos-ml": dict(type="GMRES", algebra="trilinos", preconditioner="trilinos-ml",
                        assembly="fsils", tol=1e-6, max_iter=1000, krylov=200),
}


VALVE_THRESHOLD = 0.02  # valve open: flow above this fraction of the reference peak flow
SCALARS = [  # key, label, plot marker
    ("V", "volume", "s-"), ("P", "pressure", "o-"), ("rho", "radial displacement", "D-"),
    ("SV", "stroke volume", "^-"), ("Q_in", "inlet flow", "v-"), ("Q_out", "outlet flow", "<-"),
    ("W", "stroke work", ">-"), ("EF", "ejection fraction", "x-"),
    ("dPdt+", "peak dP/dt", "+-"), ("dPdt-", "min dP/dt", "*-")]
TWIST = ("beta", "twist rate", "p-")  # helical fibers: at mid-length
FIELDS = [  # key, label, plot marker
    ("E_RR", "radial strain", "o--"), ("E_TT", "hoop strain", "s--"),
    ("E_ZZ", "axial strain", "^--"), ("J", "Jacobian, element mean", "D--"),
    ("lambda_f", "fiber stretch", "v--")]


def valve_events(t, Q, threshold):
    """Crossings of Q(t) through the threshold as (time, opening)."""
    above = Q > threshold
    return [(float(t[i] + (threshold - Q[i]) * (t[i + 1] - t[i]) / (Q[i + 1] - Q[i])),
             bool(above[i + 1])) for i in np.nonzero(above[1:] != above[:-1])[0]]


def cycle_metrics(t, V, P):
    """Stroke volume, stroke work (PV-loop area), ejection fraction and the
    extrema of dP/dt over the time points t."""
    dPdt = np.gradient(P, t)
    return {"SV": V.max() - V.min(), "W": -np.trapezoid(P, V), "EF": (V.max() - V.min()) / V.max(),
            "dPdt+": dPdt.max(), "dPdt-": dPdt.min()}


def errors(sol, ref, ref_fields, period, k, twisting=False):
    """Errors over cycle k (3D and 0D on the same time points): scalar errors,
    field errors, valve events (3D, 0D) and cycle metrics (3D, 0D). With
    twisting, also the twist rate (relative to its max over the cycle)."""
    t3 = np.asarray(sol["t"])
    idx = [int(np.argmin(np.abs(ref["t"] - t))) for t in t3]
    if not np.allclose(ref["t"][idx], t3, rtol=0, atol=1e-9):
        raise ValueError("0D and 3D time points differ")
    m = cycle_mask(t3, period, k)
    if k == 1:  # skip the start-up: the initial state (reference volume at the
        m &= t3 >= STARTUP  # initial pressure) is equilibrated in the first steps
    t = t3[m]
    on_3d = lambda key: np.interp(t3, sol["t_0d"], sol[key])[m]  # noqa: E731
    on_ref = lambda x: np.asarray(x)[idx][m]  # noqa: E731
    V, Vr = np.asarray(sol["V"])[m], on_ref(ref["V"])
    P, Pr = on_3d("P"), on_ref(ref["P"])
    rho = np.asarray(sol["rho"])[m]
    rr = np.array([np.interp(sol["R"], ref["R"], row) for row in on_ref(ref["rho"])])
    e = dict(V=np.max(np.abs(V / Vr - 1.0)), P=np.max(np.abs(P - Pr)) / (Pr.max() - Pr.min()),
             rho=np.max(np.abs(rho - rr)) / np.max(np.abs(rr)))
    if twisting:
        br = on_ref(ref["beta"])
        e["beta"] = np.max(np.abs(np.asarray(sol["beta"])[m] - br)) / np.max(np.abs(br))
    metrics = cycle_metrics(t, V, P), cycle_metrics(t, Vr, Pr)
    e.update({q: abs(metrics[0][q] / metrics[1][q] - 1.0)
              for q in ("SV", "W", "EF", "dPdt+", "dPdt-")})
    events, dt = {}, []
    for q in ("Q_in", "Q_out"):
        Q, Qr = on_3d(q), on_ref(ref[q])
        e[q] = np.max(np.abs(Q - Qr)) / np.max(np.abs(Qr))
        events[q] = [valve_events(t, x, VALVE_THRESHOLD * Qr.max()) for x in (Q, Qr)]
        if [o for _, o in events[q][0]] == [o for _, o in events[q][1]]:
            dt += [abs(a - b) for (a, _), (b, _) in zip(*events[q])]
        else:  # a valve opens or closes in one model only
            dt.append(np.inf)
    e["valves"] = 1e3 * max(dt)
    f3 = {q: np.asarray(x)[m] for q, x in sol["fields"].items()}
    fe = RS.field_errors(dict(R=sol["R"], fields=f3), {q: on_ref(x) for q, x in ref_fields.items()},
                         ref["R"], Pr, sol["order"])
    return {q: float(x) for q, x in e.items()}, fe, events, metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case")
    ap.add_argument("--nr", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--rerun", action="store_true")
    ap.add_argument("--np", type=int, default=2)
    ap.add_argument("--max-time", type=float, default=7200.0)
    ap.add_argument("--cycles", type=int, help="3D cycles to run (default: the 0D input's)")
    ap.add_argument("--compare-cycle", type=int, help="cycle compared (default: the 0D input's last)")
    ap.add_argument("--ls", choices=sorted(LINEAR_SOLVERS), default="fsils")
    ap.add_argument("--reference-ne", type=int,
                    help="compare with the direct 0D solution with this many radial elements "
                         "instead of the Richardson extrapolation (output/summary_ne<N>.json, "
                         "no figure)")
    a = ap.parse_args()
    sys.stdout.reconfigure(line_buffering=True)
    root = os.path.abspath(a.case)
    model = os.path.join(root, "chamber_cylinder_0d.json")
    with open(model) as f:
        sp = json.load(f)["simulation_parameters"]
    period, cycles = float(sp["cardiac_period"]), int(sp["number_of_cardiac_cycles"])
    k = a.compare_cycle or cycles
    out_dir = os.path.join(root, "output")
    sols = {}
    for nr in a.nr:
        d = os.path.join(root, "runs", f"struct-p2_nr{nr}")
        old = None
        if os.path.exists(os.path.join(d, "case.json")):
            with open(os.path.join(d, "case.json")) as f:
                old = json.load(f)
        if not a.run:  # post-process the existing case as it was run
            if old is None:
                print(f"3D nr = {nr}: no case")
                continue
            case = old
        else:
            case = CP.generate(model, d, nr, cycles=a.cycles,
                               linear_solver=LINEAR_SOLVERS[a.ls])
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
        if case["cycles"] < k:
            print(f"3D nr = {nr}: ran {case['cycles']} cycle(s), cycle {k} is compared")
            continue
        sols[nr] = dict(CP.post(d), order=case["order"], h=(case["Re"] - case["Ri"]) / nr)
    ref, ref_name = reference(model, out_dir, a.reference_ne)
    with open(model) as f:
        cyl = CP.find_cylinder(json.load(f))["zero_d_element_values"]
    ref_fields = RS.zero_d_fields(ref["R"], ref["rho"], ref["phi"], ref["eta"], ref["beta"], cyl,
                                  (cyl["alpha_endo"], cyl["alpha_epi"]), stresses=False)
    m = cycle_mask(ref["t"], period, k)
    print(f"{ref_name}, cycle {k}: V {ref['V'][m].min() * 1e6:.2f}-{ref['V'][m].max() * 1e6:.2f}"
          f" mL, P {ref['P'][m].min() / MMHG:.1f}-{ref['P'][m].max() / MMHG:.1f} mmHg")
    twisting = RS.has_twist(ref, cyl["length"])  # helical fibers
    scalars = SCALARS + ([TWIST] if twisting else [])
    if twisting:
        print("helical fibers: wall kinematics and strains of the 3D mid-length section "
              "(beta: local twist rate); volume, pressure and flows of the whole tube")
    summary = dict(reference=ref_name, cycle=k, twisting=twisting, three_d={})
    rows, field_rows = [], []
    for nr, sol in sols.items():
        e, fe, events, (m3, m0) = errors(sol, ref, ref_fields, period, k, twisting)
        summary["three_d"][f"nr{nr}"] = dict(
            h=sol["h"], errors=e, field_errors=fe, valve_events=events,
            metrics={q: float(x) for q, x in m3.items()})
        summary["reference_metrics"] = {q: float(x) for q, x in m0.items()}
        rows.append(RS.table_row(f"3D nr = {nr}", e, scalars + [("valves", "", "")], 10))
        field_rows.append(RS.table_row(f"3D nr = {nr}", fe, FIELDS, 10))
    print("max relative error over the cycle (valves: max |dt| [ms])\n"
          + RS.table_header(scalars + [("valves", "", "")], 10) + "\n" + "\n".join(rows))
    print("strains (relative to max |E_RR|, |E_TT|), J - 1, lambda_f - 1\n"
          + RS.table_header(FIELDS, 10) + "\n" + "\n".join(field_rows))
    os.makedirs(out_dir, exist_ok=True)
    name = f"summary_ne{a.reference_ne}.json" if a.reference_ne else "summary.json"
    with open(os.path.join(out_dir, name), "w") as f:
        json.dump(summary, f, indent=1)
    if a.reference_ne:
        return
    plot(sols, summary, ref, ref_fields, ref_name, period, k, out_dir, cyl["length"])


def plot(sols, summary, ref, ref_fields, ref_name, period, k, out_dir, L):
    """All comparisons as panels of CASE/output/comparison_cycle<k>.pdf."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    twisting = summary.get("twisting", False)
    names = ("pv", "V", "P", "Q", "dPdt", "lambda_f", "E_ES", "E_ED", "conv")
    if twisting:
        names = names[:-1] + ("twist", "conv")
        fig, axs = plt.subplots(3, 4, figsize=(20.8, 13.2))
        for x in axs.flat[len(names):]:
            x.remove()
        axs = np.array(axs.flat[:len(names)])
    else:
        fig, axs = plt.subplots(3, 3, figsize=(15.6, 13.2))
    ax = dict(zip(names, axs.flat))
    m = cycle_mask(ref["t"], period, k)
    if k == 1:  # the compared time points
        m &= ref["t"] >= STARTUP
    t0 = (k - 1) * period
    tr = ref["t"][m] - t0
    i_es = int(np.nonzero(m)[0][np.argmin(ref["V"][m])])  # end systole: minimum volume
    i_ed = int(np.nonzero(m)[0][np.argmax(ref["V"][m])])  # end diastole: maximum volume
    R = np.asarray(ref["R"]) * 1e3
    flows = (("Q_in", "-", "inlet"), ("Q_out", "--", "outlet"))
    ax["pv"].plot(ref["V"][m] * 1e6, ref["P"][m] / MMHG, "k-", lw=2.5, label=ref_name)
    ax["V"].plot(tr, ref["V"][m] * 1e6, "k-", lw=2.5, label=ref_name)
    ax["P"].plot(tr, ref["P"][m] / MMHG, "k-", lw=2.5, label=ref_name)
    for q, ls, name in flows:
        ax["Q"].plot(tr, ref[q][m] * 1e6, "k" + ls, lw=2.0, label=f"{name}, {ref_name}")
    ax["dPdt"].plot(tr, np.gradient(ref["P"], ref["t"])[m] / MMHG, "k-", lw=2.5, label=ref_name)
    if twisting:
        ax["twist"].plot(tr, np.degrees(ref["beta"][m] * L), "k-", lw=2.5, label=ref_name)
    for side, j, ls in (("endo", 0, "-"), ("epi", -1, "--")):
        ax["lambda_f"].plot(tr, ref_fields["lambda_f"][m, j], "k" + ls, lw=2.0,
                            label=f"{side}, {ref_name}")
    for key, i0, label in (("E_ES", i_es, "end systole"), ("E_ED", i_ed, "end diastole")):
        for q, ls in (("E_RR", "-"), ("E_TT", "--"), ("E_ZZ", ":")):
            ax[key].plot(R, ref_fields[q][i0], "k" + ls, lw=1.5, label=f"{q}, {ref_name}")
        ax[key].set(title=f"strain at {label} (t = {ref['t'][i0] - t0:.3f} s)")
    for i, (nr, s) in enumerate(sols.items()):
        t = np.asarray(s["t"])
        c = cycle_mask(t, period, k)
        if k == 1:
            c &= t >= STARTUP
        on = lambda key: np.interp(t, s["t_0d"], s[key])[c]  # noqa: E731
        V, P, lab, col = np.asarray(s["V"])[c], on("P"), f"3D nr{nr}", f"C{i}"
        ax["pv"].plot(V * 1e6, P / MMHG, "--", c=col, label=lab)
        ax["V"].plot(t[c] - t0, V * 1e6, "--", c=col, label=lab)
        ax["P"].plot(t[c] - t0, P / MMHG, "--", c=col, label=lab)
        for q, ls, _ in flows:
            ax["Q"].plot(t[c] - t0, on(q) * 1e6, ls, c=col, lw=1.0,
                         label=lab if q == "Q_in" else None)
        ax["dPdt"].plot(t[c] - t0, np.gradient(P, t[c]) / MMHG, "--", c=col, label=lab)
        if twisting:
            ax["twist"].plot(t[c] - t0, np.degrees(np.asarray(s["beta"])[c] * L), "--", c=col,
                             label=lab)
        lam = np.asarray(s["fields"]["lambda_f"])[c]
        for j, ls in ((0, "-"), (-1, "--")):
            ax["lambda_f"].plot(t[c] - t0, lam[:, j], ls, c=col, lw=1.0,
                                label=lab if j == 0 else None)
        for key, i0 in (("E_ES", i_es), ("E_ED", i_ed)):
            i3 = int(np.argmin(np.abs(t - ref["t"][i0])))
            for q, mk in (("E_RR", "o"), ("E_TT", "s"), ("E_ZZ", "^")):
                ax[key].plot(np.asarray(s["R"]) * 1e3, s["fields"][q][i3], mk, ms=4, mfc="none",
                             c=col, label=lab if q == "E_RR" else None)
    ax["pv"].set(xlabel="cavity volume [mL]", ylabel="cavity pressure [mmHg]", title="PV loop")
    ax["V"].set(xlabel="time in cycle [s]", ylabel="cavity volume [mL]", title="volume")
    ax["P"].set(xlabel="time in cycle [s]", ylabel="cavity pressure [mmHg]", title="pressure")
    ax["Q"].set(xlabel="time in cycle [s]", ylabel="flow [mL/s]",
                title="valve flows (inlet solid, outlet dashed)")
    ax["dPdt"].set(xlabel="time in cycle [s]", ylabel="dP/dt [mmHg/s]", title="dP/dt")
    ax["lambda_f"].set(xlabel="time in cycle [s]", ylabel="fiber stretch lambda_f",
                       title="fiber stretch (endo solid, epi dashed)")
    for key in ("E_ES", "E_ED"):
        ax[key].set(xlabel="reference radius R [mm]", ylabel="Green-Lagrange strain")
    if twisting:
        ax["twist"].set(xlabel="time in cycle [s]", ylabel="twist rate x length, beta L [deg]",
                        title="twist (3D: local rate at mid-length)")

    # mesh convergence of all errors
    a = ax["conv"]
    lev = summary["three_d"]
    if lev:
        h = np.array([r["h"] for r in lev.values()]) * 1e3
        scalars = SCALARS + ([TWIST] if twisting else [])
        for group, quantities in (("errors", scalars), ("field_errors", FIELDS)):
            for q, lab, mk in quantities:
                err = [r[group].get(q) for r in lev.values()]
                if None not in err and max(err) > 0.0:
                    a.loglog(h, err, mk, label=lab)
        e0 = next(iter(lev.values()))["errors"]["V"]
        a.loglog(h, e0 * (h / h[0]) ** 2, "k:", label="slope 2")
    a.set(xlabel="element size h = (Re - Ri)/nr [mm]", ylabel="relative error vs 0D",
          title="3D mesh convergence (struct-p2)")
    a.grid(True, which="both", alpha=0.3)
    a.legend(fontsize=6, ncol=2)
    for n, x in enumerate(axs.flat):
        if x is not a:
            x.grid(alpha=0.3)
            x.legend(fontsize=7)
        title = x.get_title()
        x.set_title("")
        x.set_title(f"({chr(97 + n)}) {title}", loc="left", fontsize=10)
    where = " (wall kinematics and strains: mid-length section)" if twisting else ""
    fig.suptitle(f"active contraction, cycle {k}: 0D ({ref_name}) vs coupled 3D "
                 f"svMultiPhysics{where}", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"comparison_cycle{k}.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
