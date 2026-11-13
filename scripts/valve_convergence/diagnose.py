"""Diagnose Newton convergence of valve blocks in svZeroDSolver test cases.

For every test case that contains a valve, run the simulation and plot, over
time: the pressures on both sides of each valve, the flow through it, its
valve status, the chamber volume (if any) and the number of Newton iterations
the generalized-alpha step needed.

Usage (from the repository root, with the compiled module on PYTHONPATH):

    python scripts/valve_convergence/diagnose.py [case ...] -o out_dir
"""

import argparse
import json
import os
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import svzerod

# Test cases that contain at least one valve block.
VALVE_CASES = [
    "valve_tanh",
    "chamber_elastance_inductor",
    "chamber_sphere",
    "chamber_sphere_straindependent_Caruel",
    "chamber_sphere_closed_loop",
    "closedLoopHeart_singleVessel_decomposed",
    "closed_loop_two_hill",
    "piecewise_Chamber_and_Valve",
]

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CASE_DIR = os.path.join(REPO, "tests", "cases")


def run_case(name, last_cycle_only=False):
    """Run one test case and return (config, dataframe, convergence dict)."""
    with open(os.path.join(CASE_DIR, name + ".json")) as f:
        cfg = json.load(f)
    sim = cfg["simulation_parameters"]
    sim["output_all_cycles"] = not last_cycle_only
    sim["output_variable_based"] = True
    sim["output_interval"] = 1
    # Keep going past non-converged steps so we can see them in the plot.
    sim["max_iter_error_to_warning"] = True

    solver = svzerod.Solver(cfg)
    solver.run()
    df = solver.get_full_result()
    conv = {
        "time": np.asarray(solver.get_step_times()),
        "iter": np.asarray(solver.get_nonlinear_iterations()),
        "residual": np.asarray(solver.get_step_residuals()),
        "equation": np.asarray(solver.get_step_residual_equations()),
    }
    return cfg, df, conv


def pivot(df):
    """Long result frame -> {variable name: array over time}, plus time axis."""
    wide = df.pivot(index="time", columns="name", values="y")
    return wide.index.to_numpy(), {c: wide[c].to_numpy() for c in wide.columns}


def valve_names(cfg):
    return [v["name"] for v in cfg.get("valves", [])]


def valve_ports(cfg, var_names):
    """For each valve, find the (p_in, p_out, q) variable names."""
    ports = {}
    for valve in cfg.get("valves", []):
        name = valve["name"]
        p_in = [v for v in var_names if v.startswith("pressure:") and v.endswith(":" + name)]
        p_out = [v for v in var_names if v.startswith("pressure:" + name + ":")]
        q = [v for v in var_names if v.startswith("flow:" + name + ":")]
        ports[name] = (
            p_in[0] if p_in else None,
            p_out[0] if p_out else None,
            q[0] if q else None,
        )
    return ports


def summarize(name, cfg, conv):
    it = conv["iter"]
    worst = it >= max(2, np.percentile(it, 95))
    blocking = Counter(conv["equation"][worst]).most_common(4)
    vtypes = sorted({v["type"] for v in cfg.get("valves", [])})
    return {
        "case": name,
        "valve_type": "+".join(vtypes),
        "n_valves": len(cfg.get("valves", [])),
        "steps": len(it),
        "mean_iter": float(it.mean()),
        "max_iter": int(it.max()),
        "total_iter": int(it.sum()),
        "n_nonconverged": int((it >= cfg["simulation_parameters"].get("maximum_nonlinear_iterations", 30)).sum()),
        "blocking_equations": blocking,
    }


def plot_case(name, cfg, df, conv, out_dir, n_cycles_shown=2):
    t_all, var = pivot(df)
    names = list(var)
    ports = valve_ports(cfg, names)
    volumes = [v for v in names if v.startswith("volume:") or v.startswith("Vc:")]
    statuses = {v: "valve_status:" + v for v in valve_names(cfg) if "valve_status:" + v in var}

    # Show only the last couple of cardiac cycles, where the solution is
    # periodic and the convergence behaviour is representative.
    period = cfg["simulation_parameters"].get("cardiac_period")
    if period is None:
        period = t_all[-1]
    t0 = max(t_all[0], t_all[-1] - n_cycles_shown * period)
    m_out = t_all >= t0
    m_it = conv["time"] >= t0

    nv = max(1, len(ports))
    n_rows = nv + (1 if volumes else 0) + 2
    fig, axes = plt.subplots(
        n_rows, 1, figsize=(9, 2.0 * n_rows), sharex=True, constrained_layout=True
    )
    axes = np.atleast_1d(axes)
    row = 0

    for valve, (p_in, p_out, q) in ports.items():
        ax = axes[row]
        if p_in:
            ax.plot(t_all[m_out], var[p_in][m_out], label="$p_\\mathrm{in}$", lw=1.2)
        if p_out:
            ax.plot(t_all[m_out], var[p_out][m_out], label="$p_\\mathrm{out}$", lw=1.2)
        ax.set_ylabel(f"{valve}\npressure")
        ax.legend(loc="upper right", fontsize=7, ncol=2)
        axq = ax.twinx()
        if q:
            axq.plot(t_all[m_out], var[q][m_out], color="0.5", lw=0.8, label="$Q$")
        axq.set_ylabel("flow", color="0.5", fontsize=8)
        axq.tick_params(axis="y", colors="0.5", labelsize=7)
        row += 1

    if statuses:
        ax = axes[row]
        for valve, col in statuses.items():
            ax.plot(t_all[m_out], var[col][m_out], lw=1.0, label=valve)
        ax.set_ylabel("valve status\n(0 open, 1 closed)")
        ax.legend(loc="upper right", fontsize=7, ncol=2)
        row += 1
    else:
        # Piecewise valves have no status variable; derive it from the sign of
        # the pressure difference, which is what the block switches on.
        ax = axes[row]
        for valve, (p_in, p_out, _) in ports.items():
            if p_in and p_out:
                ax.plot(
                    t_all[m_out],
                    (var[p_out][m_out] >= var[p_in][m_out]).astype(float),
                    lw=1.0,
                    label=valve,
                )
        ax.set_ylabel("closed\n($p_\\mathrm{out}\\geq p_\\mathrm{in}$)")
        ax.legend(loc="upper right", fontsize=7, ncol=2)
        row += 1

    if volumes:
        ax = axes[row]
        for v in volumes:
            ax.plot(t_all[m_out], var[v][m_out], lw=1.2, label=v.split(":")[-1])
        ax.set_ylabel("volume")
        ax.legend(loc="upper right", fontsize=7, ncol=2)
        row += 1

    ax = axes[row]
    ax.step(conv["time"][m_it], conv["iter"][m_it], where="post", lw=1.0, color="C3")
    ax.set_ylabel("Newton\niterations")
    ax.set_xlabel("time [s]")
    ax.set_ylim(bottom=0)
    s = summarize(name, cfg, conv)
    ax.set_title(
        f"mean {s['mean_iter']:.2f}, max {s['max_iter']}, total {s['total_iter']}",
        fontsize=8,
    )

    fig.suptitle(f"{name}  ({s['valve_type']})", fontsize=11)
    path = os.path.join(out_dir, f"convergence_{name}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path, s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="*", default=None)
    ap.add_argument("-o", "--out", default="valve_diagnostics")
    args = ap.parse_args()
    cases = args.cases or VALVE_CASES
    os.makedirs(args.out, exist_ok=True)

    summaries = []
    for name in cases:
        cfg, df, conv = run_case(name)
        path, s = plot_case(name, cfg, df, conv, args.out)
        summaries.append(s)
        print(f"{path}")

    print()
    hdr = f"{'case':42s} {'valve':16s} {'steps':>7s} {'mean':>6s} {'max':>5s} {'total':>8s}"
    print(hdr)
    print("-" * len(hdr))
    for s in summaries:
        print(
            f"{s['case']:42s} {s['valve_type']:16s} {s['steps']:7d} "
            f"{s['mean_iter']:6.2f} {s['max_iter']:5d} {s['total_iter']:8d}"
        )
    print()
    for s in summaries:
        print(f"{s['case']}: slowest equations {s['blocking_equations']}")
    with open(os.path.join(args.out, "summary.json"), "w") as f:
        json.dump(summaries, f, indent=2, default=str)


if __name__ == "__main__":
    main()
