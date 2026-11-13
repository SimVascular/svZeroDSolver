"""Robustness benchmark for the svZeroDSolver valve blocks.

Runs a physiologic left-ventricle model (see model.py) over wide ranges of
valve parameters, operating conditions and unit systems, and reports both the
cost of the non-linear solve (Newton iterations, non-converged steps) and the
accuracy of the result (stroke volume, regurgitant fraction).

Usage, from the repository root with the compiled module on PYTHONPATH:

    python scripts/valve_convergence/benchmark.py [experiment ...] -o out_dir
"""

import argparse
import contextlib
import io
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import svzerod  # noqa: E402

from model import lv_afterload, with_valve  # noqa: E402

MMHG = 1333.22  # Barye per mmHg

# Valve settings compared throughout. "tanh (tuned)" is what a careful user
# arrives at for this model in mmHg; the diode variants need no tuning.
VARIANTS = {
    "tanh": dict(valve_type="ValveTanh"),
    "piecewise": dict(valve_type="PiecewiseValve"),
    "diode": dict(valve_type="ValveDiode", eps=0.0),
    "diode-eps": dict(valve_type="ValveDiode", eps=0.1),
}


def build(variant, unit_scale=1.0, valve_numbers_fixed=False, **kwargs):
    """Build a model for one valve variant.

    ``valve_numbers_fixed`` keeps the literal numbers of the valve parameters
    when the unit system changes, which is what happens when a model is
    ported from mmHg to Barye without touching the valve block.
    """
    spec = VARIANTS[variant]
    cfg = lv_afterload(unit_scale=unit_scale, **kwargs)
    p = 1.0 if valve_numbers_fixed else unit_scale
    if spec["valve_type"] == "ValveTanh":
        params = dict(
            Rmax=kwargs.get("rmax", 1.0e4) * p,
            Rmin=kwargs.get("rmin", 0.005) * p,
            Steepness=kwargs.get("steepness", 10.0) / p,
        )
    elif spec["valve_type"] == "PiecewiseValve":
        params = dict(
            Rmax=kwargs.get("rmax", 1.0e4) * p,
            Rmin=kwargs.get("rmin", 0.005) * p,
        )
    else:
        params = dict(
            R=kwargs.get("rmin", 0.005) * p,
            eps=spec["eps"] * p,
        )
    return with_valve(cfg, spec["valve_type"], params)


def measure(cfg, unit_scale=1.0):
    """Run a configuration and return cost and accuracy metrics."""
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            solver = svzerod.Solver(cfg)
            solver.run()
            df = solver.get_full_result()
            it = np.asarray(solver.get_nonlinear_iterations())
    except Exception as exc:  # noqa: BLE001 - we want to report any failure
        return dict(ok=False, error=str(exc).split("\n")[0][:80])

    w = df.pivot(index="time", columns="name", values="y")
    t = w.index.to_numpy()
    vol = w["Vc:ventricle"].to_numpy()
    q_ao = w["flow:aortic:aorta"].to_numpy()
    q_mv = w["flow:mitral:ventricle"].to_numpy()
    p_ao = w["pressure:aortic:aorta"].to_numpy() / unit_scale

    def regurg(q):
        fwd = np.trapezoid(np.clip(q, 0.0, None), t)
        bwd = -np.trapezoid(np.clip(q, None, 0.0), t)
        return 100.0 * bwd / fwd if fwd > 0.0 else np.nan

    max_iter = cfg["simulation_parameters"].get(
        "maximum_nonlinear_iterations", 30
    )
    return dict(
        ok=True,
        mean_iter=float(it.mean()),
        max_iter=int(it.max()),
        total_iter=int(it.sum()),
        n_failed=int((it >= max_iter).sum()),
        edv=float(vol.max()),
        esv=float(vol.min()),
        sv=float(vol.max() - vol.min()),
        p_sys=float(p_ao.max()),
        p_dia=float(p_ao.min()),
        regurg_aortic=float(regurg(q_ao)),
        regurg_mitral=float(regurg(q_mv)),
    )


def row(label, r, ref_sv=None):
    if not r["ok"]:
        return f"{label:34s}  FAILED: {r['error']}"
    err = ""
    if ref_sv is not None and np.isfinite(ref_sv) and ref_sv != 0.0:
        err = f" {100 * (r['sv'] - ref_sv) / ref_sv:+7.2f}%"
    return (
        f"{label:34s} {r['mean_iter']:6.2f} {r['max_iter']:5d} "
        f"{r['total_iter']:8d} {r['n_failed']:6d} "
        f"{r['sv']:8.2f}{err} {r['regurg_aortic']:8.2f} "
        f"{r['p_dia']:5.0f}/{r['p_sys']:<5.0f}"
    )


HEADER = (
    f"{'':34s} {'mean':>6s} {'max':>5s} {'total':>8s} {'fail':>6s} "
    f"{'SV':>8s} {'SV err':>8s} {'regurg%':>8s} {'Pao':>11s}"
)


def reference_sv(**kwargs):
    """Stroke volume of the exact (non-smoothed) diode: the ground truth."""
    r = measure(build("diode", **kwargs))
    return r["sv"] if r["ok"] else np.nan


# --------------------------------------------------------------------------
# Experiments
# --------------------------------------------------------------------------


def exp_valve_parameters(out):
    """How much do the regularization parameters change cost and answer?"""
    print("\n### 1. Sensitivity to the valve regularization parameters\n")
    ref = reference_sv()
    print(f"Reference (exact diode, no regularization): SV = {ref:.2f} ml\n")
    print(HEADER)
    records = []
    for rmax in (1.0e2, 1.0e3, 1.0e4, 1.0e6, 1.0e8):
        for steep in (0.1, 1.0, 10.0, 100.0, 1000.0):
            cfg = build("tanh", rmax=rmax, steepness=steep)
            r = measure(cfg)
            print(row(f"tanh  Rmax={rmax:.0e} k={steep:<6g}", r, ref))
            records.append(dict(variant="tanh", rmax=rmax, steepness=steep, **r))
        print()
    for eps in (10.0, 1.0, 0.1, 0.01, 0.0):
        cfg = with_valve(
            lv_afterload(), "ValveDiode", dict(R=0.005, eps=eps)
        )
        r = measure(cfg)
        print(row(f"diode eps={eps:<10g}", r, ref))
        records.append(dict(variant="diode", eps=eps, **r))
    r = measure(build("piecewise"))
    print(row("piecewise", r, ref))
    records.append(dict(variant="piecewise", **r))
    return records


def exp_units(out):
    """Same physics, different unit system, valve numbers left untouched."""
    print("\n### 2. Unit system\n")
    print(
        "Each model is the identical physical problem. 'rescaled' converts "
        "every\nparameter including the valve; 'numbers kept' leaves the "
        "valve block\nuntouched, which is what porting a model from mmHg to "
        "Barye looks like.\n"
    )
    records = []
    for name, scale in (("mmHg", 1.0), ("Barye", MMHG)):
        print(f"-- {name} --")
        print(HEADER)
        ref = reference_sv(unit_scale=scale)
        for variant in ("tanh", "piecewise", "diode"):
            r = measure(build(variant, unit_scale=scale), unit_scale=scale)
            print(row(f"{variant} (rescaled)", r, ref))
            records.append(dict(units=name, variant=variant, rescaled=True, **r))
        if scale != 1.0:
            for variant in ("tanh", "piecewise", "diode"):
                r = measure(
                    build(variant, unit_scale=scale, valve_numbers_fixed=True),
                    unit_scale=scale,
                )
                print(row(f"{variant} (numbers kept)", r, ref))
                records.append(
                    dict(units=name, variant=variant, rescaled=False, **r)
                )
        print()
    return records


def exp_operating_point(out):
    """Valve parameters fixed, operating point swept: a parametric study."""
    print("\n### 3. Parametric study at fixed valve settings\n")
    sweeps = {
        "heart rate": ("period", [0.4, 0.5, 0.6, 0.8, 1.0, 1.4]),
        "contractility": ("emax", [0.5, 1.0, 2.5, 5.0, 10.0]),
        "afterload": ("rd", [0.25, 0.5, 1.0, 2.0, 4.0]),
        "preload": ("preload", [2.0, 5.0, 10.0, 20.0, 30.0]),
        "open resistance": ("rmin", [5.0e-4, 1.0e-3, 5.0e-3, 2.0e-2, 1.0e-1]),
    }
    records = []
    for title, (key, values) in sweeps.items():
        print(f"-- {title} --")
        print(HEADER)
        for value in values:
            kw = {key: value}
            ref = reference_sv(**kw)
            for variant in ("tanh", "piecewise", "diode"):
                r = measure(build(variant, **kw))
                print(row(f"{variant}  {key}={value:<8g}", r, ref))
                records.append(dict(sweep=title, param=key, value=value,
                                    variant=variant, **r))
            print()
    return records


def exp_tolerance(out):
    """Cost as a function of the absolute tolerance, in both unit systems."""
    print("\n### 4. Absolute tolerance\n")
    print(
        "The termination criterion is an absolute bound on the largest "
        "residual,\nso its meaning changes with the unit system even though "
        "the physics does not.\n"
    )
    records = []
    for name, scale in (("mmHg", 1.0), ("Barye", MMHG)):
        print(f"-- {name} --")
        print(HEADER)
        for atol in (1.0e-6, 1.0e-8, 1.0e-10):
            for variant in ("tanh", "diode"):
                r = measure(
                    build(variant, unit_scale=scale, atol=atol),
                    unit_scale=scale,
                )
                print(row(f"{variant}  atol={atol:.0e}", r))
                records.append(dict(units=name, atol=atol, variant=variant, **r))
            print()
    return records


def exp_line_search(out):
    """Does a globalized Newton iteration help?"""
    print("\n### 5. Newton line search\n")
    records = []
    print(HEADER)
    for variant in ("tanh", "diode"):
        for steps in (0, 4, 8):
            cfg = build(variant)
            cfg["simulation_parameters"]["line_search_steps"] = steps
            r = measure(cfg)
            print(row(f"{variant}  line_search_steps={steps}", r))
            records.append(dict(variant=variant, line_search=steps, **r))
        print()
    return records


def exp_time_integration(out):
    """Does the numerical damping of the generalized-alpha scheme help?"""
    print("\n### 6. Spectral radius of the time integrator\n")
    records = []
    print(HEADER)
    for variant in ("tanh", "diode"):
        for rho in (0.0, 0.25, 0.5, 0.75, 1.0):
            cfg = build(variant)
            cfg["simulation_parameters"]["rho_infty"] = rho
            r = measure(cfg)
            print(row(f"{variant}  rho_infty={rho}", r))
            records.append(dict(variant=variant, rho_infty=rho, **r))
        print()
    return records


EXPERIMENTS = {
    "valve-parameters": exp_valve_parameters,
    "units": exp_units,
    "operating-point": exp_operating_point,
    "tolerance": exp_tolerance,
    "line-search": exp_line_search,
    "time-integration": exp_time_integration,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("experiments", nargs="*", default=None)
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()
    names = args.experiments or list(EXPERIMENTS)
    all_records = {}
    for name in names:
        all_records[name] = EXPERIMENTS[name](args.out)
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, "benchmark.json"), "w") as f:
            json.dump(all_records, f, indent=2, default=str)


if __name__ == "__main__":
    main()
