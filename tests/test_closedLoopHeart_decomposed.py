"""
Test that builds a four-chamber closed-loop heart from individual blocks
(ChamberElastanceInductor, ValveTanh, BloodVessel) and compares with the
monolithic ClosedLoopHeartPulmonary block.

Matching model components:
  - Ventricular elastance: 'fourier' activation reproduces the monolithic
    25-mode Fourier series exactly
  - Pulmonary: CLHPulmonary flow-through windkessel (Q_in=Q_out) matches
    the monolithic pulmonary model
  - RA activation: 'wrapping_cosine' reproduces the monolithic cosine AA(t)
    that wraps across the cardiac cycle boundary
  - RA P-V: exponential passive mode (Kxp, Kxv, Vaso parameters) reproduces
    the monolithic nonlinear P-V relationship

All chamber and vessel parameters match the monolithic model exactly:
  - Ventricular elastance: 'fourier' activation
  - Atrial activation: 'wrapping_cosine' (cycle-boundary wrapping)
  - Atrial P-V: exponential passive mode (Kxp, Kxv, Vaso)
  - Pulmonary: CLHPulmonary flow-through windkessel

Structural limitation (~8% mean pressure/flow gap):
  In the monolithic block, each valve state directly gates the flow
  coefficient in the adjacent chamber's volume equation:
      dV_RV = valves[tricuspid] * Q_RA  -  valves[pulmonary] * Q_RV
  When the tricuspid closes (valves[tricuspid]=0), Q_RA is completely
  decoupled from dV_RV in the system matrix. With independent blocks,
  Q always contributes with coefficient ±1 regardless of the upstream
  valve state. This structural difference changes the steady-state
  operating point and cannot be eliminated by parameter matching.
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


def extract_vessel_data(result, vessel_name):
    mask = result["name"] == vessel_name
    return result[mask].sort_values("time").reset_index(drop=True)


def make_comparison_plots(mono_aorta, decomp_aorta, pdf_path):
    """Create PDF with pressure and flow comparison plots."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(
        "Monolithic ClosedLoopHeartPulmonary vs Decomposed Individual Blocks",
        fontsize=13,
    )

    T = 1.0169  # cardiac period
    mono_t = mono_aorta["time"].values % T
    decomp_t = decomp_aorta["time"].values % T

    panels = [
        (axes[0, 0], "pressure_in", "Aortic Inlet Pressure"),
        (axes[0, 1], "pressure_out", "Aortic Outlet Pressure"),
        (axes[1, 0], "flow_in", "Aortic Inlet Flow"),
        (axes[1, 1], "flow_out", "Aortic Outlet Flow"),
    ]

    for ax, col, title in panels:
        mono_vals = mono_aorta[col].values
        decomp_vals = decomp_aorta[col].values
        ax.plot(mono_t, mono_vals, "k-", lw=1.5, label="Monolithic")
        ax.plot(decomp_t, decomp_vals, "r--", lw=1.5, label="Decomposed")
        ax.set_title(title)
        ax.set_xlabel("Time in cycle [s]")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    os.makedirs(os.path.dirname(pdf_path), exist_ok=True)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"  Plot saved to {pdf_path}")


def test_decomposed_closedLoopHeart():
    """Four-chamber heart from individual blocks vs monolithic.

    Uses fourier activation for ventricles and tuned ValveTanh parameters.
    Compares aortic pressure/flow waveforms and generates PDF plots.

    Remaining differences come from:
      1. Atrial elastance model (linear vs exponential passive)
      2. Valve model (smooth tanh vs hard open/close)
      3. Pulmonary BloodVessel (capacitor stores flow) vs monolithic
         flow-through windkessel (Q_in=Q_out): In the monolithic model,
         the pulmonary RC circuit only affects the pressure relationship
         (Cp*dP_pul = Q_in - (P_pul-P_LA)/Rpd) while blood volume flows
         straight through (LA inflow = RV outflow). A standard BloodVessel
         instead has Q_in - Q_out = C*d(P_in - R*Q_in), so the capacitor
         absorbs and releases flow, introducing a phase shift between inlet
         and outlet flows. This redistributes instantaneous flow and changes
         downstream filling dynamics.
    """
    # Run both models at same time resolution
    mono_config = load_config("closedLoopHeart_singleVessel.json")
    mono_config["simulation_parameters"]["number_of_time_pts_per_cardiac_cycle"] = 1000
    mono_result = pysvzerod.simulate(mono_config)

    decomp_config = load_config("closedLoopHeart_singleVessel_decomposed.json")
    decomp_result = pysvzerod.simulate(decomp_config)

    # Extract aorta vessel data (same vessel name in both models, after Cpa)
    mono_aorta = extract_vessel_data(mono_result, "branch0_seg0")
    decomp_aorta = extract_vessel_data(decomp_result, "branch0_seg0")

    assert len(mono_aorta) > 0, "No branch0_seg0 data in monolithic result"
    assert len(decomp_aorta) > 0, "No sys_artery data in decomposed result"

    # Generate PDF comparison plots
    pdf_path = os.path.join(plot_dir, "closedLoopHeart_decomposed_comparison.pdf")
    make_comparison_plots(mono_aorta, decomp_aorta, pdf_path)

    # Compare mean hemodynamics (smooth valve gives ~5% gap from Rmax vs post_solve)
    rtol = 0.10

    print("\n" + "=" * 80)
    print("Four-chamber heart: individual blocks vs monolithic")
    print("=" * 80)

    all_pass = True
    for col, label in [
        ("pressure_in", "Aortic pressure"),
        ("flow_out", "Aortic flow"),
    ]:
        mono_mean = np.mean(mono_aorta[col].values)
        decomp_mean = np.mean(decomp_aorta[col].values)
        rel_diff = abs(mono_mean - decomp_mean) / max(abs(mono_mean), 1e-10)

        status = "PASS" if rel_diff <= rtol else "FAIL"
        if status == "FAIL":
            all_pass = False

        print(
            f"  {label:20s}: mono={mono_mean:8.2f}  decomp={decomp_mean:8.2f}"
            f"  rel_diff={rel_diff:.4f}  [{status}]"
        )

    # Pulse pressure
    mono_pp = np.ptp(mono_aorta["pressure_in"].values)
    decomp_pp = np.ptp(decomp_aorta["pressure_in"].values)
    print(f"  {'Pulse pressure':20s}: mono={mono_pp:8.2f}  decomp={decomp_pp:8.2f}")

    print()
    print("Smooth ValveTanh with Rmin/Rmax. Both at 1000 pts/cycle.")
    print("=" * 80)

    assert all_pass, (
        f"Mean hemodynamics differ by more than {rtol * 100:.0f}%. "
        "See comparison above."
    )

    # Sanity: decomposed model produces physiologically reasonable output
    mean_p = np.mean(decomp_aorta["pressure_in"].values)
    mean_q = np.mean(decomp_aorta["flow_out"].values)
    assert mean_p > 0, f"Mean aortic pressure should be positive, got {mean_p}"
    assert mean_q > 0, f"Mean aortic flow should be positive, got {mean_q}"
    assert decomp_pp > 0, f"Pulse pressure should be positive, got {decomp_pp}"
