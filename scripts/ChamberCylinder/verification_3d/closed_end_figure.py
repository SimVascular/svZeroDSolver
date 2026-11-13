"""Figure: free-length (closed) tube, 0D vs 3D along the axis.

The 0D block with fixed_length = 0 closes the tube with lids that carry the
cavity pressure (genet23 Eqs. 39-44). The 3D model closes it with a meshed end
plate (cylinder_3d.py, closed_end): an eighth of the tube with u_z = 0 on the
mid-length plane and the plate at the other end, pressure on the endocardium
and on the inner face of the plate, nothing holding the length. The plate also
holds the end ring in, so a layer near it differs from the uniform 0D state.

(a) endocardial radial displacement and (b) local axial strain du_z/dZ along
    the axis, 3D vs 0D at 120 mmHg, for tubes of 2, 4, 6 and 8 L0 (nr = 2);
(c) the same for the 6 L0 tube at nr = 1, 2, 4: the end layer is
    mesh-independent, the middle converges to the 0D block;
(d) axial strain over the pressure ramp for the three walls, 0D and 3D at
    the mid-length plane.

Usage: python closed_end_figure.py
Writes passive_isotropic_closed/output/closed_end.pdf.
"""
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
LEVELS = {1: "#86b6ef", 2: "#3987e5", 4: "#1c5cab", 8: "#0d366b"}  # blue ramp: nr 1, 2, 4, 8
LENGTHS = [("passive_isotropic_closed_L2", 2, "#f6b48f"), ("passive_isotropic_closed_L4", 4, "#eb6834"),
           ("passive_isotropic_closed", 6, "#b8461a"), ("passive_isotropic_closed_L8", 8, "#6e2a0f")]
WALLS = [("passive_isotropic_closed", "isotropic", "#2a78d6", "o"),
         ("passive_anisotropic_closed", "one fiber family", "#eb6834", "s"),
         ("passive_anisotropic_full_closed", "full Holzapfel-Ogden", "#1baf7a", "d")]


def reference(case):
    """0D Richardson reference (ne, 2 ne) and the 0D input values."""
    root = os.path.join(HERE, case)
    m0 = C3.read_0d(os.path.join(root, "chamber_cylinder_0d.json"))
    ne = int(m0["values"]["num_elements"])
    out = os.path.join(root, "output")
    ref = RS.richardson(RS.run_0d(m0["config"], m0["name"], ne, out),
                        RS.run_0d(m0["config"], m0["name"], 2 * ne, out))
    return ref, m0


def complete(d):
    """The run folder d holds all load steps."""
    if not os.path.exists(os.path.join(d, "case.json")):
        return False
    with open(os.path.join(d, "case.json")) as f:
        n = json.load(f)["n_steps"]
    return len(glob.glob(os.path.join(d, "results", "result_*.vtu"))) == n


def axial_profiles(case, nr):
    """Distance from the plate (in wall thicknesses H), endocardial rho
    (3D/0D - 1) and local axial strain (3D/0D - 1) at the last step."""
    d = os.path.join(HERE, case, "runs", f"struct-p2_nr{nr}")
    if not complete(d):
        return None
    with open(os.path.join(d, "case.json")) as f:
        c = json.load(f)
    m = C3.mesh_closed_cylinder(c["Ri"], c["Re"], c["mesh_length"], c["nr"], c["plate"], c["order"],
                                c["ntheta"], c["nz"], c["nz_plate"], c["n_in"])
    X = m["points"]
    files = sorted(glob.glob(os.path.join(d, "results", "result_*.vtu")),
                   key=lambda fn: int(os.path.basename(fn)[7:-4]))
    _, pdata = C3._read_vtu(files[-1])
    u = pdata["Displacement"][:, :3]
    ref, m0 = reference(case)
    T, _ = C3.pressure_ramp(m0["t"], m0["P"])
    i = int(np.argmin(np.abs(ref["t"] - T)))
    rho0, eps0 = float(ref["rho"][i][0]), float(ref["eps"][i])
    Ri, Re, Lm = c["Ri"], c["Re"], c["mesh_length"]
    R = np.hypot(X[:, 0], X[:, 1])
    wall = (R > Ri - 1e-9) & (X[:, 2] < Lm + 1e-9)
    endo = wall & (np.abs(R - Ri) < 1e-9)
    Z = np.unique(np.round(X[wall, 2], 12))
    s = c["nz"] * c["order"]
    Z = Z[np.rint(Z / Lm * s).astype(int) % c["order"] == 0]  # vertex levels (mid-edge levels zigzag)
    rho = np.array([(np.hypot(*(X[sel, :2] + u[sel, :2]).T) - Ri).mean()
                    for sel in (endo & (np.abs(X[:, 2] - z) < 1e-12) for z in Z)])
    uz = np.array([u[wall & (np.abs(X[:, 2] - z) < 1e-12), 2].mean() for z in Z])
    H = Re - Ri
    return dict(d=(Lm - Z) / H, rho=rho / rho0 - 1.0, eps=np.gradient(uz, Z) / eps0 - 1.0,
                L=Lm / H)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "legend.frameon": False})
    fig, axs = plt.subplots(2, 2, figsize=(11.0, 8.0))
    style = dict(lw=2.0)

    def floor(a):
        a.set_yscale("log")
        a.set_ylim(1e-5, 2.0)
        a.set_xlabel("distance from the end plate / wall thickness $H$")
        a.grid(True, which="major", color="0.92", lw=0.6)

    # (a), (b): tube length, nr = 2
    for case, n, c in LENGTHS:
        p = axial_profiles(case, 2)
        if p is None:
            continue
        lab = f"${n}L_0$" + (" (used)" if n == 6 else "")
        axs[0, 0].plot(p["d"], np.abs(p["rho"]), color=c, label=lab, **style)
        axs[0, 1].plot(p["d"], np.abs(p["eps"]), color=c, label=lab, **style)
        for a, k in ((axs[0, 0], "rho"), (axs[0, 1], "eps")):
            a.plot(p["d"][0], abs(p[k][0]), "o", color=c, ms=8, mfc="white", mew=2.0)
    axs[0, 0].set(title=r"(a) isotropic, nr = 2: endocardial $\rho$, $|$3D/0D $- 1|$",
                  ylabel=r"$|\rho_\mathrm{3D}/\rho_\mathrm{0D} - 1|$ at 120 mmHg")
    axs[0, 1].set(title=r"(b) isotropic, nr = 2: local axial strain, $|$3D/0D $- 1|$",
                  ylabel=r"$|\partial u_z/\partial Z\,/\,\varepsilon_\mathrm{0D} - 1|$ at 120 mmHg")
    # (c): mesh levels, 6 L0
    for nr in (1, 2, 4):
        p = axial_profiles("passive_isotropic_closed", nr)
        if p is None:
            continue
        axs[1, 0].plot(p["d"], np.abs(p["rho"]), color=LEVELS[nr], label=f"nr = {nr}", **style)
        axs[1, 0].plot(p["d"][0], abs(p["rho"][0]), "o", color=LEVELS[nr], ms=8, mfc="white", mew=2.0)
    axs[1, 0].set(title=r"(c) isotropic, $6L_0$: endocardial $\rho$ per mesh level",
                  ylabel=r"$|\rho_\mathrm{3D}/\rho_\mathrm{0D} - 1|$ at 120 mmHg")
    for a in (axs[0, 0], axs[0, 1], axs[1, 0]):
        floor(a)
        a.plot([], [], "o", color=MUTED, ms=8, mfc="white", mew=2.0, label="mid-length plane")
        a.legend(fontsize=8, loc="upper right")
    # (d): axial strain over the ramp, three walls
    a = axs[1, 1]
    mmHg = 133.322
    for case, label, c, m in WALLS:
        s_path = os.path.join(HERE, case, "output", "summary.json")
        if not os.path.exists(s_path):
            continue
        with open(s_path) as f:
            s = json.load(f)
        ref = s["reference"]
        a.plot(np.asarray(ref["P"]) / mmHg, np.asarray(ref["eps"]) * 100, color=c, lw=2.0,
               label=f"{label}, 0D")
        nr = max(int(k[2:]) for k in s["three_d"])
        r = s["three_d"][f"nr{nr}"]
        a.plot(np.asarray(r["P"]) / mmHg, np.asarray(r["eps"]) * 100, m, color=c, ms=7,
               mfc="none", mew=1.5, label=f"{label}, 3D nr = {nr}")
    a.set(title="(d) axial strain at the mid-length plane over the ramp",
          xlabel="cavity pressure [mmHg]", ylabel=r"axial strain $\varepsilon$ [%]")
    a.grid(True, which="major", color="0.92", lw=0.6)
    a.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out = os.path.join(HERE, "passive_isotropic_closed", "output", "closed_end.pdf")
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
