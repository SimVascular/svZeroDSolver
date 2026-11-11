"""Figure: 0D vs 3D for a fixed-length cylinder with helical fibers.

The 3D model is as free as the 0D: the full cylinder, u_z = 0 on both ends
(fixed length), twist and radial expansion free (rigid-body modes pinned, see
cylinder_3d.py). No 0D data enters the 3D model.

(a) helix angle through the wall: mirroring theta or Z flips its sign, so the
    eighth model's symmetry planes do not apply;
(b) Cauchy stresses of a cross-section of the 0D solution at the end of the
    ramp: the 0D sections carry sigma_thz, the free 3D ends cannot;
(c, d) radial displacement of the endocardium and rotation of the mid-wall
    along the axis at the end of the ramp, 0D and 3D per mesh level;
(e, f) 0D end shear and end twist against the helix angle.

Usage: python helical_figure.py CASE_DIR
Writes CASE_DIR/output/helical_discrepancy.pdf.
"""
import argparse
import copy
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cylinder_3d as C3  # noqa: E402
import run_study as RS  # noqa: E402

INK, MUTED = "#0b0b0b", "#8a8984"
SHEAR = "#eb6834"  # categorical slot 2
LEVELS = ["#a9c9f2", "#6aa2e8", "#2a78d6", "#1a4f93"]  # one hue, light to dark: nr 1, 2, 4, 8
FAMILY = {"±α (endo α, epi −α)": "#1baf7a", "uniform α": "#6250d6"}  # slots 3, 7


def material_values(values):
    """The 0D material parameters as check_materials.stress_0d takes them."""
    vals = {k: float(values.get(k, 0.0)) for k in ("C1", "C2", "C3", "C4", "C5", "C6")}
    vals.update({k: float(x) for k, x in values.items()
                 if k.startswith("ho_") or k == "bulk_modulus"})
    return vals


def section_stresses(sol, i, values, n=65):
    """Cauchy stresses (deformed cylindrical components) of the 0D solution
    sol at time index i on n radii."""
    from scipy.interpolate import CubicSpline
    from check_materials import stress_0d
    Ri, Re = float(values["Ri"]), float(values["Re"])
    a0, a1 = float(values.get("alpha_endo", 0.0)), float(values.get("alpha_epi", 0.0))
    R = np.linspace(Ri, Re, n)
    sp = {k: CubicSpline(sol["R"], sol[k][i]) for k in ("rho", "phi", "eta")}
    vals = material_values(values)
    beta = float(sol["beta"][i])
    out = {k: np.zeros(n) for k in ("rr", "tt", "zz", "tz")}
    for j, Rj in enumerate(R):
        r = Rj + sp["rho"](Rj)
        F = np.array([[1.0 + sp["rho"](Rj, 1), 0.0, 0.0],
                      [r * sp["phi"](Rj, 1), r / Rj, r * beta],
                      [sp["eta"](Rj, 1), 0.0, 1.0]])
        a = np.radians(a0 + (a1 - a0) * (Rj - Ri) / (Re - Ri))
        f = np.array([0.0, np.cos(a), np.sin(a)])
        sig = F @ stress_0d(F, f, np.array([1.0, 0.0, 0.0]), vals) @ F.T / np.linalg.det(F)
        for k, (p, q) in (("rr", (0, 0)), ("tt", (1, 1)), ("zz", (2, 2)), ("tz", (1, 2))):
            out[k][j] = sig[p, q]
    return R, out


def axial_profiles(case_dir):
    """Radial displacement of the endocardium and rotation of the mid-wall
    (theta-averaged) against Z at the last step of a 3D case."""
    with open(os.path.join(case_dir, "case.json")) as f:
        case = json.load(f)
    mesh = C3.mesh_cylinder(case["Ri"], case["Re"], case["mesh_length"], case["nr"],
                            case["sector"], case["order"], case["ntheta"], case["nz"])
    X = mesh["points"]
    files = sorted(glob.glob(os.path.join(case_dir, "results", "result_*.vtu")),
                   key=lambda fn: int(os.path.basename(fn)[7:-4]))
    _, pdata = C3._read_vtu(files[-1])
    x = X + pdata["Displacement"][:, :3]
    s = case["nr"] * case["order"]
    lev = np.rint((np.hypot(X[:, 0], X[:, 1]) - case["Ri"]) / (case["Re"] - case["Ri"]) * s).astype(int)
    kz = np.rint(X[:, 2] / case["mesh_length"] * case["nz"] * case["order"]).astype(int)
    rho = np.hypot(x[:, 0], x[:, 1]) - np.hypot(X[:, 0], X[:, 1])
    dth = np.arctan2(x[:, 1], x[:, 0]) - np.arctan2(X[:, 1], X[:, 0])
    dth = (dth + np.pi) % (2.0 * np.pi) - np.pi
    Z = np.unique(kz) * case["mesh_length"] / (case["nz"] * case["order"])
    endo, mid = lev == 0, lev == s // 2
    prof = lambda m, f: np.array([f[m & (kz == k)].mean() for k in np.unique(kz)])  # noqa: E731
    return Z, prof(endo, rho), prof(mid, dth)


def finished_levels(root, suffix="", formulation="struct-p2"):
    """Mesh levels nr with complete 3D results in root/runs/<formulation>_nr<nr><suffix>."""
    out = []
    for d in sorted(glob.glob(os.path.join(root, "runs", f"{formulation}_nr*{suffix}"))):
        tail = os.path.basename(d)[len(formulation) + 3:]
        tail = tail[:len(tail) - len(suffix)] if suffix else tail
        if not tail.isdigit() or not os.path.exists(os.path.join(d, "case.json")):
            continue
        with open(os.path.join(d, "case.json")) as f:
            n_steps = json.load(f)["n_steps"]
        if len(glob.glob(os.path.join(d, "results", "result_*.vtu"))) == n_steps:
            out.append((int(tail), d))
    return sorted(out)


def sweep(config, name, out_dir, angles, T, ne=16):
    """0D end shear max|sigma_thz| / max sigma_thth and end twist beta L (deg)
    at the end of the ramp for symmetric and uniform helix angles."""
    res = {}
    for fam in FAMILY:
        shear, twist = [], []
        for a in angles:
            cfg = copy.deepcopy(config)
            v = next(b for b in cfg["vessels"] if b["vessel_name"] == name)["zero_d_element_values"]
            a = float(a)
            v["alpha_endo"], v["alpha_epi"] = (a, -a) if fam.startswith("±") else (a, a)
            sol = RS.run_0d(cfg, name, ne, out_dir)
            i = int(np.argmin(np.abs(sol["t"] - T)))
            _, sig = section_stresses(sol, i, v)
            shear.append(np.max(np.abs(sig["tz"])) / np.max(np.abs(sig["tt"])))
            twist.append(np.degrees(sol["beta"][i] * v["length"]))
        res[fam] = (np.array(shear), np.array(twist))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", help="folder with chamber_cylinder_0d.json and runs/struct-p2_nr*")
    ap.add_argument("--suffix", default="", help="run folder suffix, e.g. _damp0.001 (run_study --end-damping)")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    root = os.path.abspath(a.case)
    m0 = C3.read_0d(os.path.join(root, "chamber_cylinder_0d.json"))
    v, out_dir = m0["values"], os.path.join(root, "output")
    T, Pmax = C3.pressure_ramp(m0["t"], m0["P"])
    ne = int(v.get("num_elements", 10))
    ref = RS.richardson(RS.run_0d(m0["config"], m0["name"], ne, out_dir),
                        RS.run_0d(m0["config"], m0["name"], 2 * ne, out_dir))
    i = int(np.argmin(np.abs(ref["t"] - T)))
    L, Ri, Re = v["length"], v["Ri"], v["Re"]

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "legend.frameon": False, "lines.linewidth": 2})
    fig, ax = plt.subplots(2, 3, figsize=(11, 6.8), constrained_layout=True)
    mm = 1e3

    # (a) helix angle
    R = np.linspace(Ri, Re, 65)
    alpha = v["alpha_endo"] + (v["alpha_epi"] - v["alpha_endo"]) * (R - Ri) / (Re - Ri)
    ax[0, 0].plot(R * mm, alpha, color=INK)
    ax[0, 0].axhline(0.0, color=MUTED, lw=0.8)
    ax[0, 0].set(xlabel="R [mm]", ylabel="helix angle α [deg]",
                 title="(a) fiber helix angle through the wall")
    ax[0, 0].text(0.04, 0.06, "mirroring θ or Z maps α → −α:\nno symmetry planes, full cylinder",
                  transform=ax[0, 0].transAxes, color=MUTED, fontsize=8)

    # (b) section stresses of the 0D solution
    Rs, sig = section_stresses(ref, i, v)
    for k, lab, c, ls in (("tt", "σθθ", INK, "-"), ("zz", "σzz", MUTED, "-"),
                          ("rr", "σrr", MUTED, ":"), ("tz", "σθz", SHEAR, "-")):
        ax[0, 1].plot(Rs * mm, sig[k] / 1e3, color=c, ls=ls, label=lab)
    ax[0, 1].axhline(0.0, color=MUTED, lw=0.8)
    lo, hi = ax[0, 1].get_ylim()
    ax[0, 1].set_ylim(lo, hi + 0.45 * (hi - lo))  # room for the legend
    ax[0, 1].set(xlabel="R [mm]", ylabel="Cauchy stress [kPa]",
                 title=f"(b) 0D cross-section at {Pmax / C3.MMHG:.0f} mmHg")
    ax[0, 1].legend(loc="upper left", ncol=2)
    ax[0, 1].text(0.97, 0.80, "every 0D section carries σθz;\na free 3D end cannot",
                  transform=ax[0, 1].transAxes, color=MUTED, fontsize=8, ha="right", va="top")

    # (c, d) along the axis, 0D and 3D per mesh level
    rho_endo = float(ref["rho"][i][0])
    ax[0, 2].axhline(0.0, color=INK, ls="--", label="0D")
    # rotations relative to mid-length: the rigid rotation is free in both models
    ax[1, 0].plot([0.0, L * mm], np.degrees(ref["beta"][i] * L) * np.array([-0.5, 0.5]), color=INK,
                  ls="--", label="0D")
    for k, (nr, d) in enumerate(finished_levels(root, a.suffix)):
        Z, rho, th = axial_profiles(d)
        c = LEVELS[min(k, len(LEVELS) - 1)]
        ax[0, 2].plot(Z * mm, 100 * (rho / rho_endo - 1.0), color=c, marker="o", ms=3, label=f"3D nr = {nr}")
        ax[1, 0].plot(Z * mm, np.degrees(th - np.interp(0.5 * L, Z, th)), color=c, marker="o", ms=3,
                      label=f"3D nr = {nr}")
    for axx in (ax[0, 2], ax[1, 0]):
        axx.legend(loc="best")
        axx.set_xlabel("Z [mm]")
    ax[0, 2].set(ylabel="endo radial displacement: 3D / 0D − 1 [%]",
                 title="(c) radial expansion along the axis")
    ax[1, 0].set(ylabel="mid-wall rotation, relative to mid-length [deg]", title="(d) twist along the axis")

    # (e, f) 0D against the helix angle
    angles = np.arange(0, 91, 15)
    res = sweep(m0["config"], m0["name"], out_dir, angles, T)
    for fam, c in FAMILY.items():
        ax[1, 1].plot(angles, 100 * res[fam][0], color=c, marker="o", ms=4, label=fam)
        ax[1, 2].plot(angles, res[fam][1], color=c, marker="o", ms=4, label=fam)
    ax[1, 1].set(xlabel="helix angle α [deg]", ylabel="max |σθz| / max |σθθ| [%]",
                 title="(e) 0D section shear vs helix angle")
    ax[1, 2].set(xlabel="helix angle α [deg]", ylabel="end twist βL [deg]",
                 title="(f) 0D end twist vs helix angle")
    ax[1, 2].axhline(0.0, color=MUTED, lw=0.8)
    ax[1, 1].set_ylim(0.0, 1.4 * ax[1, 1].get_ylim()[1])  # room for the legend
    ax[1, 1].legend(loc="upper right")
    ax[1, 2].legend(loc="lower right")
    fig.suptitle("Fixed-length cylinder with helical fibers: 0D vs 3D at "
                 f"{Pmax / C3.MMHG:.0f} mmHg\n3D: full cylinder, both ends u_z = 0 "
                 "(as the 0D), twist and radial expansion free", fontsize=10)
    path = os.path.join(out_dir, "helical_discrepancy.pdf")
    fig.savefig(path)
    print(path)


if __name__ == "__main__":
    main()
