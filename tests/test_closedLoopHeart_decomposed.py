"""
Test that builds a four-chamber closed-loop heart from individual blocks
and compares with the monolithic ClosedLoopHeartPulmonary block.
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import pysvzerod

this_file_dir = os.path.abspath(os.path.dirname(__file__))
cases_dir = os.path.join(this_file_dir, "cases")
plot_dir = os.path.join(this_file_dir, "plots")

def load_config(filename):
    with open(os.path.join(cases_dir, filename)) as f:
        return json.load(f)


def get_var(result, name):
    """Extract a variable-based output as sorted arrays (time, y)."""
    d = result[result["name"] == name].sort_values("time")
    return d["time"].values, d["y"].values


def get_vessel(result, vessel_name):
    d = result[result["name"] == vessel_name].sort_values("time")
    return d


def test_decomposed_closedLoopHeart():
    """Four-chamber heart from individual blocks vs monolithic."""

    # ---- run both models at same resolution with variable-based output ----
    shared = {"number_of_time_pts_per_cardiac_cycle": 1000,
              "number_of_cardiac_cycles": 4,
              "output_all_cycles": False}

    mono_config = load_config("closedLoopHeart_singleVessel.json")
    mono_config["simulation_parameters"].update(shared)
    mono_config["simulation_parameters"]["output_variable_based"] = True
    mono_r = pysvzerod.simulate(mono_config)

    decomp_config = load_config("closedLoopHeart_singleVessel_decomposed.json")
    decomp_config["simulation_parameters"].update(shared)
    decomp_config["simulation_parameters"]["output_variable_based"] = True
    decomp_r = pysvzerod.simulate(decomp_config)

    # also run without variable-based for vessel-level comparison
    mono_config2 = load_config("closedLoopHeart_singleVessel.json")
    mono_config2["simulation_parameters"].update(shared)
    mono_v = pysvzerod.simulate(mono_config2)

    decomp_config2 = load_config("closedLoopHeart_singleVessel_decomposed.json")
    decomp_config2["simulation_parameters"].update(shared)
    decomp_config2["simulation_parameters"]["output_variable_based"] = False
    decomp_v = pysvzerod.simulate(decomp_config2)

    T = 1.0169

    # ---- mapping between monolithic and decomposed variable names ----
    chamber_map = [
        ("RA", "V_RA:CLH", "Vc:RA"),
        ("RV", "V_RV:CLH", "Vc:RV"),
        ("LA", "V_LA:CLH", "Vc:LA"),
        ("LV", "V_LV:CLH", "Vc:LV"),
    ]
    pressure_map = [
        ("P_RA", "pressure:J_heart_inlet:CLH", "pressure:J_return:RA"),
        ("P_RV", "P_RV:CLH", "pressure:TV:RV"),
        ("P_LA", "P_LA:CLH", "pressure:J_pul:LA"),
        ("P_LV", "P_LV:CLH", "pressure:MV:LV"),
    ]
    flow_map = [
        ("Q_RA (tricuspid)", "Q_RA:CLH", "flow:RA:TV"),
        ("Q_RV (pulmonary)", "Q_RV:CLH", "flow:RV:PV"),
        ("Q_LA (mitral)",    "Q_LA:CLH", "flow:LA:MV"),
        ("Q_LV (aortic)",   "Q_LV:CLH", "flow:LV:AV"),
    ]

    os.makedirs(plot_dir, exist_ok=True)

    def save_fig(fig, name):
        path = os.path.join(plot_dir, name)
        fig.savefig(path)
        plt.close(fig)
        print(f"  Saved {path}")

    mono_aorta = get_vessel(mono_v, "branch0_seg0")
    decomp_aorta = get_vessel(decomp_v, "branch0_seg0")

    # ---- 1: Aortic pressure and flow ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Aortic Pressure and Flow", fontsize=14)
    for ax, col, title in [
        (axes[0, 0], "pressure_in", "Aortic Inlet Pressure"),
        (axes[0, 1], "pressure_out", "Aortic Outlet Pressure"),
        (axes[1, 0], "flow_in", "Aortic Inlet Flow"),
        (axes[1, 1], "flow_out", "Aortic Outlet Flow"),
    ]:
        ax.plot(mono_aorta["time"].values , mono_aorta[col].values,
                "k-", lw=1.5, label="Monolithic")
        ax.plot(decomp_aorta["time"].values , decomp_aorta[col].values,
                "r--", lw=1.5, label="Decomposed")
        ax.set_title(title); ax.set_xlabel("Time [s]")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "decomp_01_aortic.pdf")

    # ---- 2: PV loops for all 4 chambers ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Pressure-Volume Loops", fontsize=14)
    for ax, (label, m_vol, d_vol) in zip(axes.flat, chamber_map):
        _, m_v = get_var(mono_r, m_vol)
        _, d_v = get_var(decomp_r, d_vol)
        m_p_name = [p for p in pressure_map if p[0].endswith(label)][0]
        _, m_p = get_var(mono_r, m_p_name[1])
        _, d_p = get_var(decomp_r, m_p_name[2])
        n = min(len(m_v), len(m_p), len(d_v), len(d_p))
        ax.plot(m_v[:n], m_p[:n], "k-", lw=1.5, label="Monolithic")
        ax.plot(d_v[:n], d_p[:n], "r--", lw=1.5, label="Decomposed")
        ax.set_title(f"{label} PV Loop"); ax.set_xlabel("Volume"); ax.set_ylabel("Pressure")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "decomp_02_pv_loops.pdf")

    # ---- 3: Valve flows ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Valve Flows", fontsize=14)
    for ax, (label, m_name, d_name) in zip(axes.flat, flow_map):
        m_t, m_q = get_var(mono_r, m_name)
        d_t, d_q = get_var(decomp_r, d_name)
        ax.plot(m_t , m_q, "k-", lw=1.5, label="Monolithic")
        ax.plot(d_t , d_q, "r--", lw=1.5, label="Decomposed")
        ax.set_title(label); ax.set_xlabel("Time [s]"); ax.set_ylabel("Flow")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "decomp_03_valve_flows.pdf")

    # ---- 4: Valve dQ/dt (flow rate of change) ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Valve dQ/dt (flow rate of change)", fontsize=14)
    for ax, (label, m_name, d_name) in zip(axes.flat, flow_map):
        m_t, m_q = get_var(mono_r, m_name)
        d_t, d_q = get_var(decomp_r, d_name)
        dt_m = np.diff(m_t); dt_m[dt_m == 0] = 1e-10
        dt_d = np.diff(d_t); dt_d[dt_d == 0] = 1e-10
        m_dq = np.diff(m_q) / dt_m
        d_dq = np.diff(d_q) / dt_d
        ax.plot(m_t[:-1] , m_dq, "k-", lw=1, label="Monolithic", alpha=0.8)
        ax.plot(d_t[:-1] , d_dq, "r--", lw=1, label="Decomposed", alpha=0.8)
        ax.set_title(label); ax.set_xlabel("Time [s]"); ax.set_ylabel("dQ/dt")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "decomp_05_valve_dQdt.pdf")

    # ---- 5: Chamber volumes ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Chamber Volumes", fontsize=14)
    for ax, (label, m_name, d_name) in zip(axes.flat, chamber_map):
        m_t, m_v = get_var(mono_r, m_name)
        d_t, d_v = get_var(decomp_r, d_name)
        ax.plot(m_t , m_v, "k-", lw=1.5, label="Monolithic")
        ax.plot(d_t , d_v, "r--", lw=1.5, label="Decomposed")
        ax.set_title(f"{label} Volume"); ax.set_xlabel("Time [s]"); ax.set_ylabel("Volume")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "decomp_06_volumes.pdf")

    # ---- numerical comparison ----
    mono_aorta = get_vessel(mono_v, "branch0_seg0")
    decomp_aorta = get_vessel(decomp_v, "branch0_seg0")

    rtol = 0.10
    print("\n" + "=" * 80)
    print("Four-chamber heart: individual blocks vs monolithic (1 cycle, 1000 pts)")
    print("=" * 80)

    all_pass = True
    for col, label in [("pressure_in", "Aortic pressure"), ("flow_out", "Aortic flow")]:
        mm = np.mean(mono_aorta[col].values)
        dm = np.mean(decomp_aorta[col].values)
        rd = abs(mm - dm) / max(abs(mm), 1e-10)
        status = "PASS" if rd <= rtol else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  {label:20s}: mono={mm:8.2f}  decomp={dm:8.2f}  rel_diff={rd:.4f}  [{status}]")

    mono_pp = np.ptp(mono_aorta["pressure_in"].values)
    decomp_pp = np.ptp(decomp_aorta["pressure_in"].values)
    print(f"  {'Pulse pressure':20s}: mono={mono_pp:8.2f}  decomp={decomp_pp:8.2f}")
    print("=" * 80)

    assert all_pass
    assert np.mean(decomp_aorta["pressure_in"].values) > 0
    assert np.mean(decomp_aorta["flow_out"].values) > 0

test_decomposed_closedLoopHeart()