"""Sketch: boundary conditions of the 0D and the 3D helical-fiber cylinder.

Top row: longitudinal section (Z right, R up) of the wall in the 0D model and
in the 3D model (full cylinder, both ends u_z = 0, as free as the 0D), and
why they differ. Middle row: the 0D section shear through the wall and the 3D
twist and radial expansion along the axis (finest finished 3D mesh). Bottom:
options.

Usage: python helical_sketch.py CASE_DIR
Writes CASE_DIR/output/helical_bcs_sketch.pdf.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cylinder_3d as C3  # noqa: E402
import run_study as RS  # noqa: E402
from helical_figure import INK, LEVELS, MUTED, SHEAR, axial_profiles, finished_levels, section_stresses  # noqa: E402

WALL, BAND = "#e4e3df", "#b5b4ae"
Z0, Z1, R0, R1 = 0.0, 8.0, 1.6, 4.4  # sketch coordinates of the wall
FS = 8


def wall(ax, title, caption):
    """Wall section with the cavity pressure, the axis, a title and a caption."""
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((Z0, R0), Z1 - Z0, R1 - R0, facecolor=WALL, edgecolor=BAND, lw=1))
    ax.text(0.5 * (Z0 + Z1), R1 + 0.12, "epicardium", color=MUTED, fontsize=FS - 1, ha="center", va="bottom")
    ax.plot([Z0 - 1.0, Z1 + 1.0], [-0.4, -0.4], color=MUTED, lw=0.8, ls="-.")
    ax.text(Z0 - 1.0, -0.3, "axis", color=MUTED, fontsize=FS - 1, va="bottom")
    for z in np.linspace(0.8, Z1 - 0.8, 6):  # cavity pressure on the endocardium
        ax.annotate("", xy=(z, R0 - 0.05), xytext=(z, R0 - 0.8),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1))
    ax.text(0.5 * (Z0 + Z1), R0 - 1.15, "cavity pressure P on the endocardium", color=MUTED,
            fontsize=FS - 1, ha="center", va="center")
    ax.text(Z0 - 1.6, 7.3, title, fontsize=FS + 2, fontweight="bold", va="top")
    ax.text(Z0 - 1.6, -1.0, caption, fontsize=FS, va="top", color=INK, linespacing=1.4)
    ax.set_xlim(Z0 - 1.7, Z1 + 3.0)
    ax.set_ylim(-4.2, 7.4)
    ax.set_aspect("equal")
    ax.axis("off")


def shear_symbols(ax, z, tz, color):
    """sigma_thz across the wall: dot = out of the page, cross = into it,
    size by magnitude."""
    R = np.linspace(R0 + 0.3, R1 - 0.3, len(tz))
    s = np.abs(tz) / np.max(np.abs(tz))
    for r, si, t in zip(R, s, tz):
        if si < 0.1:
            continue
        ax.scatter([z], [r], s=110 * si, facecolors="white", edgecolors=color, lw=1.2, zorder=4)
        ax.scatter([z], [r], s=(14 if t > 0 else 30) * si, marker="o" if t > 0 else "x",
                   color=color, lw=1.2, zorder=5)


def slip_end(ax, z, side):
    """Rollers: u_z = 0, free in the plane of the end."""
    from matplotlib.patches import Circle
    for r in np.linspace(R0 + 0.3, R1 - 0.3, 5):
        ax.add_patch(Circle((z + side * 0.18, r), 0.18, facecolor="white", edgecolor=INK, lw=1))
    ax.plot([z + side * 0.4] * 2, [R0 - 0.15, R1 + 0.15], color=INK, lw=1.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("--suffix", default="", help="run folder suffix, e.g. _damp0.001 (run_study --end-damping)")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    root = os.path.abspath(a.case)
    m0 = C3.read_0d(os.path.join(root, "chamber_cylinder_0d.json"))
    v, out_dir = m0["values"], os.path.join(root, "output")
    T, Pmax = C3.pressure_ramp(m0["t"], m0["P"])
    ne = int(v.get("num_elements", 10))
    ref = RS.richardson(RS.run_0d(m0["config"], m0["name"], ne, out_dir),
                        RS.run_0d(m0["config"], m0["name"], 2 * ne, out_dir))
    i = int(np.argmin(np.abs(ref["t"] - T)))
    _, sig = section_stresses(ref, i, v, n=7)
    Rf, sigf = section_stresses(ref, i, v)
    L, mm = v["length"], 1e3
    levels = finished_levels(root, a.suffix)

    plt.rcParams.update({"font.size": FS, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                         "legend.frameon": False, "lines.linewidth": 2})
    fig = plt.figure(figsize=(12, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(3, 3, height_ratios=[1.45, 0.9, 0.42])
    sk = [fig.add_subplot(gs[0, k]) for k in range(3)]
    pl = [fig.add_subplot(gs[1, k]) for k in range(3)]
    tx = fig.add_subplot(gs[2, :])
    tx.axis("off")

    # --- 0D: Z-homogeneous, every cross-section carries sigma_thz(R)
    ax = sk[0]
    wall(ax, "0D model",
         "r = R + ρ(R),  θ = Θ + βZ,  z = Z: no ends, the\n"
         "same state at every Z (an infinitely long tube).\n"
         "Every cross-section carries the shear σθz(R)\n"
         "(● out of, ✕ into the page); zero net torque sets β.")
    for z in (2.0, 6.0):
        ax.plot([z, z], [R0, R1], color=INK, lw=1, ls="--")
        shear_symbols(ax, z, sig["tz"], INK)
    ax.text(6.35, 0.5 * (R0 + R1), "σθz(R)", fontsize=FS, va="center", color=INK)

    # --- 3D: full cylinder, both ends u_z = 0, otherwise free
    ax = sk[1]
    wall(ax, "3D model",
         "Full cylinder, both ends u_z = 0 (fixed length, as\n"
         "the 0D), twist and radial expansion free; in-plane\n"
         "dashpots on both ends hold the rigid-body modes (no\n"
         "force at rest). The free ends carry no σθz: end layers.")
    for z0 in (Z0, Z1 - 2.4):
        ax.add_patch(Rectangle((z0, R0), 2.4, R1 - R0, facecolor=SHEAR, alpha=0.2, lw=0))
    ax.text(Z0 + 1.2, 0.5 * (R0 + R1), "end\nlayer", fontsize=FS, ha="center", va="center", color=SHEAR)
    ax.text(Z1 - 1.2, 0.5 * (R0 + R1), "end\nlayer", fontsize=FS, ha="center", va="center", color=SHEAR)
    slip_end(ax, Z0, -1)
    slip_end(ax, Z1, 1)
    ax.text(Z1 + 0.6, 0.5 * (R0 + R1), "u_z = 0\nσθz = 0", fontsize=FS, va="center", color=INK)

    # --- why
    ax = sk[2]
    ax.axis("off")
    ax.text(0.0, 0.93, "Why they differ", fontsize=FS + 2, fontweight="bold", va="top",
            transform=ax.transAxes)
    ax.text(0.0, 0.80,
            "Helical fibers load every cross-section with\n"
            "the azimuthal shear σθz(R): positive near the\n"
            "endocardium, negative near the epicardium,\n"
            "zero net torque (here up to 18% of σθθ).\n\n"
            "The 0D has no ends, so this costs nothing. A\n"
            "3D end that only keeps u_z = 0 has σθz = 0, so\n"
            "near each end the wall twists and expands\n"
            "differently. This tube is only about four wall\n"
            "thicknesses long: for the small ±60° twist the\n"
            "end layers reach the middle.\n\n"
            "This is a modeling difference: it does not shrink\n"
            "with mesh refinement, only with a longer tube.",
            fontsize=FS + 0.5, va="top", transform=ax.transAxes, linespacing=1.4)

    # --- middle row: 0D section shear, 3D twist and radial expansion along Z
    ax = pl[0]
    ax.plot(Rf * mm, sigf["tz"] / 1e3, color=SHEAR, label="σθz")
    ax.plot(Rf * mm, sigf["tt"] / 1e3, color=MUTED, lw=1.2, label="σθθ (for scale)")
    ax.axhline(0.0, color=MUTED, lw=0.8)
    ax.set(xlabel="R [mm]", ylabel="0D Cauchy stress [kPa]",
           title=f"0D section stress, {Pmax / C3.MMHG:.0f} mmHg, ±60° fibers")
    ax.legend(loc="center right")
    rho_endo = float(ref["rho"][i][0])
    # rotations relative to mid-length: the rigid rotation is free in both models
    pl[1].plot([0.0, L * mm], np.degrees(ref["beta"][i] * L) * np.array([-0.5, 0.5]), color=INK,
               ls="--", label="0D")
    pl[2].axhline(0.0, color=INK, ls="--", label="0D")
    if levels:
        nr, d = levels[-1]
        Z, rho, th = axial_profiles(d)
        c = LEVELS[min(len(levels) - 1, len(LEVELS) - 1)]
        pl[1].plot(Z * mm, np.degrees(th - np.interp(0.5 * L, Z, th)), color=c, marker="o", ms=3,
                   label=f"3D nr = {nr}")
        pl[2].plot(Z * mm, 100 * (rho / rho_endo - 1.0), color=c, marker="o", ms=3, label=f"3D nr = {nr}")
    pl[1].set(xlabel="Z [mm]", ylabel="mid-wall rotation, relative to mid-length [deg]",
              title="twist along the axis")
    pl[2].set(xlabel="Z [mm]", ylabel="endo radial displacement: 3D / 0D − 1 [%]",
              title="radial expansion along the axis")
    for ax in pl[1:]:
        ax.legend(loc="best")

    tx.text(0.0, 1.0, "Options without 0D data in the 3D model", fontsize=FS + 2, fontweight="bold",
            va="top")
    tx.text(0.0, 0.72,
            "1. Compare away from the ends (mid-length section) and report the end layers as a "
            "model difference.\n"
            "2. A longer 3D tube (L / wall thickness ↑): the mid-length state approaches the 0D, "
            "at the cost of larger meshes.\n"
            "3. Circumferential or axial fibers only (σθz = 0): free ends and symmetry planes are "
            "exact (the passive and active cases so far).",
            fontsize=FS + 0.5, va="top", linespacing=1.6)
    path = os.path.join(out_dir, "helical_bcs_sketch.pdf")
    fig.savefig(path)
    print(path)


if __name__ == "__main__":
    main()
