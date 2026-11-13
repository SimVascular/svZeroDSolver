"""Figures comparing the valve formulations on the benchmark LV model.

Produces two figures:

``valve_comparison.png``
    One cardiac cycle of the benchmark model solved with each valve block,
    showing ventricular and aortic pressure, ventricular volume, the valve
    state and the Newton iterations each time step needed.

``valve_robustness.png``
    Maps of solver cost and stroke-volume error over the regularization
    parameters of ValveTanh, next to the single parameter of ValveDiode.

Usage, from the repository root with the compiled module on PYTHONPATH:

    python scripts/valve_convergence/compare.py -o out_dir
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import svzerod  # noqa: E402

from benchmark import build, measure  # noqa: E402
from model import lv_afterload, with_valve  # noqa: E402

VARIANTS = [
    ("ValveTanh", dict(Rmax=1.0e4, Rmin=0.005, Steepness=10.0), "C0", "-", 1.2),
    ("PiecewiseValve", dict(Rmax=1.0e4, Rmin=0.005), "C1", "-", 3.0),
    ("ValveDiode", dict(R=0.005, eps=0.0), "C2", "--", 1.4),
]


def run(valve_type, params):
    cfg = with_valve(lv_afterload(), valve_type, params)
    solver = svzerod.Solver(cfg)
    solver.run()
    w = solver.get_full_result().pivot(index="time", columns="name", values="y")
    return (
        w,
        np.asarray(solver.get_step_times()),
        np.asarray(solver.get_nonlinear_iterations()),
    )


def figure_cycle(out_dir):
    fig, axes = plt.subplots(
        4, 1, figsize=(9, 8.5), sharex=True, constrained_layout=True
    )
    for valve_type, params, color, ls, lw in VARIANTS:
        w, t_it, it = run(valve_type, params)
        t = w.index.to_numpy()
        label = valve_type + (
            f" (eps={params['eps']:g})" if "eps" in params else ""
        )
        axes[0].plot(t, w["pressure:ventricle:aortic"], color=color, lw=lw,
                     ls=ls, label=label)
        axes[0].plot(t, w["pressure:aortic:aorta"], color=color, lw=lw,
                     ls=":")
        axes[1].plot(t, w["Vc:ventricle"], color=color, lw=lw, ls=ls,
                     label=label)
        status = (
            w[f"valve_status:aortic"]
            if "valve_status:aortic" in w
            else (w["pressure:aortic:aorta"] >= w["pressure:ventricle:aortic"]
                  ).astype(float)
        )
        axes[2].plot(t, status, color=color, lw=lw, ls=ls, label=label)
        axes[3].step(t_it, it, where="post", color=color, lw=lw, ls=ls,
                     label=f"{label}: {it.sum()} total, max {it.max()}")

    axes[0].set_ylabel("pressure [mmHg]")
    axes[0].set_title(
        "ventricle (solid/dashed) and aorta (dotted); PiecewiseValve and "
        "ValveDiode coincide", fontsize=9, loc="left"
    )
    axes[0].legend(fontsize=8)
    axes[1].set_ylabel("LV volume [ml]")
    axes[2].set_ylabel("aortic valve status\n(0 open, 1 closed)")
    axes[3].set_ylabel("Newton iterations")
    axes[3].set_xlabel("time [s]")
    axes[3].legend(fontsize=8)
    axes[3].set_ylim(bottom=0)
    period = 0.8
    axes[3].set_xlim(t[-1] - period, t[-1])
    fig.suptitle("Benchmark LV model, last cardiac cycle", fontsize=12)
    path = os.path.join(out_dir, "valve_comparison.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def figure_robustness(out_dir):
    steeps = np.logspace(-1, 3, 9)
    rmaxes = np.logspace(2, 8, 7)
    cost = np.full((len(rmaxes), len(steeps)), np.nan)
    err = np.full_like(cost, np.nan)
    ref = measure(build("diode"))["sv"]
    for i, rmax in enumerate(rmaxes):
        for j, steep in enumerate(steeps):
            r = measure(build("tanh", rmax=rmax, steepness=steep))
            if r["ok"]:
                cost[i, j] = r["total_iter"]
                err[i, j] = 100.0 * abs(r["sv"] - ref) / ref

    epss = np.logspace(-3, 1, 9)
    d_cost, d_err = [], []
    for eps in epss:
        r = measure(
            with_valve(lv_afterload(), "ValveDiode", dict(R=0.005, eps=eps))
        )
        d_cost.append(r["total_iter"])
        d_err.append(100.0 * abs(r["sv"] - ref) / ref)
    r0 = measure(build("diode"))

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    ext = [np.log10(steeps[0]), np.log10(steeps[-1]),
           np.log10(rmaxes[0]), np.log10(rmaxes[-1])]
    for ax, data, title, cmap in (
        (axes[0, 0], cost, "ValveTanh: total Newton iterations", "viridis"),
        (axes[1, 0], np.log10(np.maximum(err, 1e-3)),
         "ValveTanh: |stroke volume error| [log10 %]", "magma"),
    ):
        im = ax.imshow(data, origin="lower", extent=ext, aspect="auto",
                       cmap=cmap)
        fig.colorbar(im, ax=ax)
        ax.set_xlabel("log10 Steepness")
        ax.set_ylabel("log10 Rmax")
        ax.set_title(title, fontsize=10)

    axes[0, 1].semilogx(epss, d_cost, "o-", color="C2", label="ValveDiode")
    axes[0, 1].axhline(r0["total_iter"], color="C2", ls=":",
                       label=f"eps = 0: {r0['total_iter']}")
    axes[0, 1].axhline(np.nanmin(cost), color="C0", ls="--",
                       label=f"best ValveTanh: {int(np.nanmin(cost))}")
    axes[0, 1].set_xlabel("eps [mmHg]")
    axes[0, 1].set_ylabel("total Newton iterations")
    axes[0, 1].set_title("ValveDiode: cost vs. regularization", fontsize=10)
    axes[0, 1].legend(fontsize=8)

    axes[1, 1].loglog(epss, np.maximum(d_err, 1e-4), "o-", color="C2")
    axes[1, 1].set_xlabel("eps [mmHg]")
    axes[1, 1].set_ylabel("|stroke volume error| [%]")
    axes[1, 1].set_title(
        "ValveDiode: error vanishes as eps -> 0", fontsize=10
    )
    axes[1, 1].grid(True, which="both", alpha=0.3)

    fig.suptitle(
        "Cost and accuracy over the valve regularization parameters",
        fontsize=12,
    )
    path = os.path.join(out_dir, "valve_robustness.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="valve_diagnostics")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    print(figure_cycle(args.out))
    print(figure_robustness(args.out))


if __name__ == "__main__":
    main()
