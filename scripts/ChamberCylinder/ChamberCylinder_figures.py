# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Reproduce the Genet (2023) Fig. 5-9 comparisons for the ChamberCylinder block.

Each figure reads the *digitized* Genet traces stored next to this script
(``genet23_fig{5..9}_digitized.csv``), runs the ChamberCylinder model for
the matching case(s), and plots both in the paper's own layout: cavity pressure,
cavity volume, PV loop and ventricular torsion, with the twin kPa/mmHg pressure
axis. Solid = digitized Genet; dashed = our model. The digitized CSVs (not the
copyrighted figures) are the only paper data kept in the repo.

  Fig 5  baseline case (single trace)
  Fig 6  aspect ratio   L/2Ri  = 1.0 / 1.5 / 2.0
  Fig 7  internal volume       = 58 / 65 / 72 mL
  Fig 8  wall volume           = 117 / 130 / 143 mL
  Fig 9  fiber angle  +-alpha  = 30 / 45 / 60 / 75 / 90 deg

Run (from the repository root, after building the Python module):
  PYTHONPATH=build/python python scripts/ChamberCylinder/ChamberCylinder_figures.py
Writes fig5.png ... fig9.png to scripts/ChamberCylinder/output/.
"""
import os
import sys
import csv
import collections

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import svzerod

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ChamberCylinder_genet_exact as G  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
kPa = 1e-3
MMHG_PER_KPA = 7.50062          # 14 kPa -> 105 mmHg (paper's right-hand axis)
GAMMA_TABLE1 = 70.0             # genet23 Table 1 viscosity [Pa s]

# Spectral colors sampled from the paper renders; variations use them in order.
CMAP = {'red': (213/255, 62/255, 80/255), 'orange': (248/255, 157/255, 91/255),
        'yellow': (241/255, 235/255, 146/255), 'green': (155/255, 209/255, 162/255),
        'blue': (46/255, 136/255, 189/255)}


# ne = 6 is within 0.2 mL and 0.05 deg of ne = 12; 12 cycles are periodic to
# 0.01 mL.
def run(overrides=None, ne=6, ncycle=12, gamma=GAMMA_TABLE1):
    """Run the genet23 baseline with optional vessel-parameter overrides
    (geometry or fiber angle) and wall viscosity gamma; return the last cycle."""
    cfg = G.build(ne=ne, ncycle=ncycle, P_vs=G.P_VS)
    vv = cfg["vessels"][1]["zero_d_element_values"]
    vv["gamma"] = gamma
    if overrides:
        vv.update(overrides)
        if {"Ri", "length"} & set(overrides):
            cfg["initial_condition"]["volume:ventricle"] = (
                np.pi * vv["Ri"] ** 2 * vv["length"])
    cfg["simulation_parameters"]["absolute_tolerance"] = 1e-6
    try:
        out = svzerod.simulate(cfg)
    except Exception as e:  # report and skip a case that fails
        print("  FAIL", str(e)[:70]); return None

    def s(n):
        d = out[out["name"] == n].sort_values("time")
        return d["time"].to_numpy(), d["y"].to_numpy()

    t, V = s("volume:ventricle")
    _, Pv = s("pressure:ventricle:aortic")
    _, b = s("beta:ventricle")
    L = vv["length"]
    m = t >= t.max() - 0.8 - 1e-9
    return dict(t=(t[m] - t[m].min()) * 1000, V=V[m] * 1e6, P=Pv[m] * kPa,
                tw=b[m] * L * 180 / np.pi)


def geom(aspect, Vint, Vwall):
    """Cylinder geometry from aspect ratio L/2Ri, internal and wall volume."""
    Ri = (Vint / (2 * np.pi * aspect)) ** (1 / 3)
    L = 2 * Ri * aspect
    return dict(Ri=Ri, Re=np.sqrt(Ri ** 2 + Vwall / (np.pi * L)), length=L)


# Baseline internal and wall volumes (65 and 130 mL, genet23 Table 1).
_d = G.build()["vessels"][1]["zero_d_element_values"]
RI0, RE0, L0 = _d["Ri"], _d["Re"], _d["length"]
VINT0 = np.pi * RI0 ** 2 * L0
VWALL0 = np.pi * L0 * (RE0 ** 2 - RI0 ** 2)

FIGURES = {
    'fig6': dict(vlabel="aspect ratio\n$L/2R_i$", vol=(0, 140), pv=(70, 140), cases=[
        ('1.0', 'red', geom(1.0, VINT0, VWALL0)),
        ('1.5', 'yellow', geom(1.5, VINT0, VWALL0)),
        ('2.0', 'blue', geom(2.0, VINT0, VWALL0))]),
    'fig7': dict(vlabel="internal\nvolume [mL]", vol=(0, 160), pv=(50, 160), cases=[
        ('58', 'red', geom(1.5, 58e-6, VWALL0)),
        ('65', 'yellow', geom(1.5, 65e-6, VWALL0)),
        ('72', 'blue', geom(1.5, 72e-6, VWALL0))]),
    'fig8': dict(vlabel="wall\nvolume [mL]", vol=(0, 140), pv=(70, 140), cases=[
        ('117', 'red', geom(1.5, VINT0, 117e-6)),
        ('130', 'yellow', geom(1.5, VINT0, 130e-6)),
        ('143', 'blue', geom(1.5, VINT0, 143e-6))]),
    'fig9': dict(vlabel="fiber angle\n$\\pm\\alpha$ [deg]", vol=(0, 140), pv=(70, 140), cases=[
        ('30', 'red', dict(alpha_endo=30, alpha_epi=-30)),
        ('45', 'orange', dict(alpha_endo=45, alpha_epi=-45)),
        ('60', 'yellow', dict(alpha_endo=60, alpha_epi=-60)),
        ('75', 'green', dict(alpha_endo=75, alpha_epi=-75)),
        ('90', 'blue', dict(alpha_endo=90, alpha_epi=-90))]),
}


# --------------------------------------------------------------------------- #
#  digitized-trace loading + plotting helpers
# --------------------------------------------------------------------------- #
def load_long(key):
    """Load a long-format digitized CSV: dict[(panel, variation)] = (x, y)."""
    tmp = collections.defaultdict(list)
    with open(os.path.join(HERE, f"genet23_{key}_digitized.csv")) as f:
        for r in csv.reader(f):
            if not r or r[0].startswith("#") or r[0] == "panel":
                continue
            tmp[(r[0], r[1])].append((float(r[2]), float(r[3])))
    out = {}
    for k, pts in tmp.items():
        pts.sort()
        out[k] = (np.array([p[0] for p in pts]), np.array([p[1] for p in pts]))
    return out


def gapbreak(x, y, gap=16.0):
    """Insert NaN where the digitized trace has a gap, so a partly-occluded curve
    is drawn only where it was actually visible (no straight-line bridges)."""
    xn, yn = [x[0]], [y[0]]
    for i in range(1, len(x)):
        if x[i] - x[i - 1] > gap:
            xn.append(np.nan); yn.append(np.nan)
        xn.append(x[i]); yn.append(y[i])
    return np.array(xn), np.array(yn)


def _far(t, s, gap):
    j = np.searchsorted(s, t); d = np.full(len(t), 1e9)
    for k, ti in enumerate(t):
        c = []
        if j[k] < len(s): c.append(abs(s[j[k]] - ti))
        if j[k] > 0: c.append(abs(ti - s[j[k] - 1]))
        if c: d[k] = min(c)
    return d > gap


def paper_pvloop(D, var, gap=16.0):
    """Reconstruct the paper PV loop parametrically from digitized V(t), P(t)."""
    tv, V = D[('volume', var)]; tp, P = D[('pressure', var)]
    t = np.linspace(0, 800, 700)
    Vi = np.interp(t, tv, V); Pi = np.interp(t, tp, P)
    bad = _far(t, tv, gap) | _far(t, tp, gap)
    Vi[bad] = np.nan; Pi[bad] = np.nan
    return Vi, Pi


def _style_time(a, ylim, ylabel, yticks):
    a.set_xlim(0, 800); a.set_ylim(*ylim)
    a.set_xticks(np.arange(0, 801, 100)); a.set_yticks(yticks)
    a.set_xlabel("time [ms]"); a.set_ylabel(ylabel)
    a.grid(True, axis="x", ls=":", lw=0.6, color="0.7")


def _add_mmhg(a):
    r = a.twinx(); lo, hi = a.get_ylim()
    r.set_ylim(lo * MMHG_PER_KPA, hi * MMHG_PER_KPA)
    r.set_ylabel("cavity pressure [mmHg]")


def _four_panels():
    fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.6))
    return fig, ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]


def _finish_axes(P, Vp, PV, Tw, F):
    ptk = np.arange(0, 15, 2); ttk = np.arange(-5, 31, 5)
    _style_time(P, (0, 14), "cavity pressure [kPa]", ptk); _add_mmhg(P)
    _style_time(Vp, F['vol'], "cavity volume [mL]",
                np.arange(0, F['vol'][1] + 1, 20))
    _style_time(Tw, (-5, 30), "ventricular torsion [$\\degree$]", ttk)
    PV.set_xlim(*F['pv']); PV.set_ylim(0, 14); PV.set_yticks(ptk)
    PV.set_xlabel("cavity volume [mL]"); PV.set_ylabel("cavity pressure [kPa]")
    PV.grid(True, ls=":", lw=0.6, color="0.85"); _add_mmhg(PV)


def overlay(key, gamma=GAMMA_TABLE1, tag="", note=""):
    F = FIGURES[key]; D = load_long(key)
    print(f"{key}{tag}: running {len(F['cases'])} case(s)")
    res = [run(ov, gamma=gamma) for _, _, ov in F['cases']]
    fig, P, Vp, PV, Tw = _four_panels()
    for (lab, cname, _), r in zip(F['cases'], res):
        col = CMAP[cname]
        for a, panel in [(P, 'pressure'), (Vp, 'volume'), (Tw, 'torsion')]:
            x, y = gapbreak(*D[(panel, lab)]); a.plot(x, y, '-', color=col, lw=1.9)
        Vi, Pi = paper_pvloop(D, lab); PV.plot(Vi, Pi, '-', color=col, lw=1.9)
        if r:
            P.plot(r['t'], r['P'], '--', color=col, lw=1.4)
            Vp.plot(r['t'], r['V'], '--', color=col, lw=1.4)
            Tw.plot(r['t'], r['tw'], '--', color=col, lw=1.4)
            PV.plot(r['V'], r['P'], '--', color=col, lw=1.4)
            print(f"  {lab:>4}: EDV={r['V'].max():.0f} ESV={r['V'].min():.0f} "
                  f"Ppk={r['P'].max():.1f} tw_pk={r['tw'].max():.0f}")
    _finish_axes(P, Vp, PV, Tw, F)
    handles = [Line2D([0], [0], color=CMAP[c], lw=3) for _, c, _ in F['cases']]
    handles.append(Line2D([0], [0], color="0.2", lw=1.7, ls="--"))
    labels = [lab for lab, _, _ in F['cases']] + ["ours"]
    fig.legend(handles, labels, title=F['vlabel'], loc="center left",
               bbox_to_anchor=(1.0, 0.5), frameon=True, fontsize=11, title_fontsize=11)
    fig.suptitle(f"Genet (2023) {key.upper()} — solid: digitized paper   "
                 f"dashed: ChamberCylinder{note}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, f"{key}{tag}.png"), dpi=120, bbox_inches="tight")
    plt.close(); print(f"  saved {key}{tag}.png")


def fig5(gamma=GAMMA_TABLE1, tag="", note=""):
    """Baseline case vs the digitized Fig 5 (single trace, black)."""
    import io
    path = os.path.join(HERE, "genet23_fig5_digitized.csv")
    with open(path) as f:
        body = "".join(ln for ln in f if not ln.lstrip().startswith("#"))
    d = np.genfromtxt(io.StringIO(body), delimiter=",", names=True)
    r = run(gamma=gamma); print(f"fig5{tag}: EDV={r['V'].max():.0f} ESV={r['V'].min():.0f} "
                                f"Ppk={r['P'].max():.1f} tw_pk={r['tw'].max():.0f}")
    fig, P, Vp, PV, Tw = _four_panels()
    ours = CMAP['red']
    P.plot(d["time_ms"], d["pressure_kPa"], 'k-', lw=1.9)
    P.plot(r['t'], r['P'], '--', color=ours, lw=1.5)
    Vp.plot(d["time_ms"], d["volume_mL"], 'k-', lw=1.9)
    Vp.plot(r['t'], r['V'], '--', color=ours, lw=1.5)
    PV.plot(d["volume_mL"], d["pressure_kPa"], 'k-', lw=1.9)
    PV.plot(r['V'], r['P'], '--', color=ours, lw=1.5)
    Tw.plot(d["time_ms"], d["torsion_deg"], 'k-', lw=1.9)
    Tw.plot(r['t'], r['tw'], '--', color=ours, lw=1.5)
    _finish_axes(P, Vp, PV, Tw, dict(vol=(0, 140), pv=(70, 140)))
    handles = [Line2D([0], [0], color='k', lw=3),
               Line2D([0], [0], color=ours, lw=1.7, ls="--")]
    fig.legend(handles, ["Genet 2023", "ours"], loc="center left",
               bbox_to_anchor=(1.0, 0.5), frameon=True, fontsize=11)
    fig.suptitle(f"Genet (2023) FIG5 — solid: digitized paper   dashed: ChamberCylinder{note}",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, f"fig5{tag}.png"), dpi=120, bbox_inches="tight")
    plt.close(); print(f"  saved fig5{tag}.png")


if __name__ == "__main__":
    # Full genet23 Table-1 damping (gamma = mu = 70 Pa s).
    fig5()
    for k in ("fig6", "fig7", "fig8", "fig9"):
        overlay(k)
