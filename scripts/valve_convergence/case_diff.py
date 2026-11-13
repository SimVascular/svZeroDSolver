"""Compare a test case against an earlier revision of its configuration.

Runs the case as it is now and as it was at a git revision, and reports both
the cost of the non-linear solve and how far the solution moved. Used to
document what switching the test suite from `ValveTanh` and `PiecewiseValve`
to `ValveDiode` changed.

Usage, from the repository root with the compiled module on PYTHONPATH:

    python scripts/valve_convergence/case_diff.py [case ...] --rev HEAD
"""

import argparse
import contextlib
import io
import json
import os
import subprocess
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CASE_DIR = os.path.join("tests", "cases")

import svzerod  # noqa: E402


def valve_cases():
    """Every test case that has at least one valve block."""
    out = []
    for name in sorted(os.listdir(os.path.join(REPO, CASE_DIR))):
        if not name.endswith(".json"):
            continue
        with open(os.path.join(REPO, CASE_DIR, name)) as f:
            try:
                cfg = json.load(f)
            except json.JSONDecodeError:
                continue
        if cfg.get("valves"):
            out.append(name[: -len(".json")])
    return out


def config_at(case, rev):
    """The configuration of `case` at a git revision, or None if absent."""
    path = f"{CASE_DIR}/{case}.json"
    r = subprocess.run(
        ["git", "show", f"{rev}:{path}"], cwd=REPO, capture_output=True, text=True
    )
    if r.returncode != 0:
        return None
    return json.loads(r.stdout)


def run(cfg):
    """Run a configuration, returning the solution and the solver cost."""
    cfg = json.loads(json.dumps(cfg))
    cfg["simulation_parameters"]["output_variable_based"] = True
    cfg["simulation_parameters"]["max_iter_error_to_warning"] = True
    max_iter = cfg["simulation_parameters"].get(
        "maximum_nonlinear_iterations", 30
    )
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        solver = svzerod.Solver(cfg)
        solver.run()
        df = solver.get_full_result()
        it = np.asarray(solver.get_nonlinear_iterations())
    wide = df.pivot(index="time", columns="name", values="y")
    return wide, dict(
        mean=float(it.mean()),
        max=int(it.max()),
        total=int(it.sum()),
        failed=int((it >= max_iter).sum()),
    )


def compare(old, new):
    """How far each variable moved, separating transient from sustained.

    Every measure is relative to the variable's own amplitude over the
    reported cycle, so a spike at the instant a valve switches is not
    confused with a shift of the whole trace.
    """
    shared = [c for c in new.columns if c in old.columns]
    moves = {}
    for name in shared:
        a, b = old[name].to_numpy(), new[name].to_numpy()
        if a.shape != b.shape:
            continue
        scale = max(np.abs(a).max(), 1e-30)
        span = max(a.max() - a.min(), 1e-30)
        moves[name] = dict(
            peak=float(np.abs(b - a).max() / scale),
            mean=float(abs(b.mean() - a.mean()) / scale),
            span=float(abs((b.max() - b.min()) - (a.max() - a.min())) / span),
            frac=float((np.abs(b - a) / scale > 0.01).mean()),
        )
    only_new = [c for c in new.columns if c not in old.columns]
    only_old = [c for c in old.columns if c not in new.columns]
    return moves, only_new, only_old


def kind(name):
    return name.split(":", 1)[0]


def integrated(old, new):
    """Chamber volumes before and after, as end-diastolic and stroke volume."""
    out = []
    for name in new.columns:
        if kind(name) not in ("Vc", "volume") or name not in old.columns:
            continue
        a, b = old[name].to_numpy(), new[name].to_numpy()
        out.append((name, a.max(), a.max() - a.min(), b.max(), b.max() - b.min()))
    return out


def switch_timing(old, new):
    """How many output samples the valve open/close instants moved by."""
    out = []
    for name in new.columns:
        if kind(name) != "valve_status" or name not in old.columns:
            continue
        a = old[name].to_numpy() > 0.5
        b = new[name].to_numpy() > 0.5
        out.append((name, int((a != b).sum()), len(a)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="*", default=None)
    ap.add_argument("--rev", default="HEAD", help="revision to compare against")
    ap.add_argument("--top", type=int, default=3, help="variables listed per case")
    args = ap.parse_args()
    cases = args.cases or valve_cases()

    hdr = (
        f"{'case':42s} {'mean':>11s} {'max':>9s} {'total':>15s} {'failed':>11s}"
    )
    print("PERFORMANCE   (before -> after)")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for case in cases:
        old_cfg = config_at(case, args.rev)
        with open(os.path.join(REPO, CASE_DIR, case + ".json")) as f:
            new_cfg = json.load(f)
        if old_cfg is None:
            print(f"{case:42s} {'(new case)':>11s}")
            continue
        old_sol, old_cost = run(old_cfg)
        new_sol, new_cost = run(new_cfg)
        print(
            f"{case:42s} "
            f"{old_cost['mean']:5.2f}->{new_cost['mean']:<5.2f} "
            f"{old_cost['max']:4d}->{new_cost['max']:<4d} "
            f"{old_cost['total']:7d}->{new_cost['total']:<7d} "
            f"{old_cost['failed']:5d}->{new_cost['failed']:<5d}"
        )
        rows.append((case, old_sol, new_sol, old_cost, new_cost))

    print()
    print("SOLUTION   (relative to each variable's own amplitude)")
    hdr2 = (
        f"{'case':42s} {'peak':>9s} {'cycle mean':>11s} {'amplitude':>10s} "
        f"{'samples >1%':>12s}"
    )
    print(hdr2)
    print("-" * len(hdr2))
    detail = []
    case_solutions = {}
    for case, old_sol, new_sol, _, _ in rows:
        moves, only_new, only_old = compare(old_sol, new_sol)
        if not moves:
            continue
        physics = {k: v for k, v in moves.items() if kind(k) != "valve_status"}
        if not physics:
            continue
        worst = {
            k: max(m[k] for m in physics.values())
            for k in ("peak", "mean", "span", "frac")
        }
        print(
            f"{case:42s} {worst['peak']:9.1%} {worst['mean']:11.1%} "
            f"{worst['span']:10.1%} {worst['frac']:12.1%}"
        )
        detail.append((case, moves, only_new, only_old))
        case_solutions[case] = (old_sol, new_sol)

    print()
    print("WORST VARIABLES PER CASE")
    for case, moves, only_new, only_old in detail:
        print(f"\n{case}")
        physics = {k: v for k, v in moves.items() if kind(k) != "valve_status"}
        worst = sorted(physics.items(), key=lambda kv: -kv[1]["peak"])[: args.top]
        for name, m in worst:
            print(
                f"      peak {m['peak']:7.1%}  mean {m['mean']:7.1%}  "
                f"amplitude {m['span']:7.1%}  >1% on {m['frac']:5.1%} of samples"
                f"   {name}"
            )
        for name, moved, total in switch_timing(*case_solutions[case]):
            print(
                f"      valve state differs on {moved} of {total} output "
                f"samples   {name}"
            )
        if only_new or only_old:
            print(f"      variables added: {only_new}  removed: {only_old}")

    print()
    print("CHAMBER VOLUMES   (before -> after)")
    for case, old_sol, new_sol, _, _ in rows:
        vols = integrated(old_sol, new_sol)
        if not vols:
            continue
        print(f"\n{case}")
        for name, edv_a, sv_a, edv_b, sv_b in vols:
            d_edv = 100 * (edv_b - edv_a) / abs(edv_a) if edv_a else float("nan")
            d_sv = 100 * (sv_b - sv_a) / abs(sv_a) if sv_a else float("nan")
            print(
                f"      {name:28s} EDV {edv_a:12.6g} -> {edv_b:12.6g} "
                f"({d_edv:+6.2f}%)   SV {sv_a:12.6g} -> {sv_b:12.6g} "
                f"({d_sv:+6.2f}%)"
            )


if __name__ == "__main__":
    main()
