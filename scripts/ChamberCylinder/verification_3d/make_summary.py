"""Summary of the ChamberCylinder 0D vs 3D verification as a PDF.

Collects the results of all studies (CASE/output/summary.json of run_study.py
and run_active_study.py, zero_d_convergence.json, and the helical-fiber
analyses computed here from the run folders) into tables, includes the
comparison figures, and compiles summary/verification_summary.tex with
pdflatex. Missing results show as "pending".

Usage: python make_summary.py [--no-compile]
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "summary")
PENDING = r"\textit{pending}"


def load(path):
    try:
        with open(os.path.join(HERE, path)) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def e(x, digits=1):
    """Number in LaTeX scientific notation."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "--"
    if x == 0:
        return "0"
    m, p = f"{x:.{digits}e}".split("e")
    return rf"${m}\cdot10^{{{int(p)}}}$"


def order(a, b, ha, hb):
    """Observed convergence order between two levels."""
    if a and b and a > 0 and b > 0:
        return math.log(a / b) / math.log(ha / hb)
    return None


def table(header, rows, caption, label, align=None):
    align = align or "l" + "r" * (len(header) - 1)
    lines = [r"\begin{table}[H]\centering\small", rf"\caption{{{caption}}}\label{{{label}}}",
             (r"\resizebox{\linewidth}{!}{" if len(header) > 8 else "")
             + rf"\begin{{tabular}}{{{align}}}", r"\toprule", " & ".join(header) + r" \\", r"\midrule"]
    lines += [" & ".join(str(c) for c in r) + r" \\" for r in rows]
    lines += [r"\bottomrule", r"\end{tabular}" + ("}" if len(header) > 8 else ""), r"\end{table}"]
    return "\n".join(lines)


def level_rows(summary, keys, group="errors", extra=()):
    """Rows nr, cells, errors (with the observed order in parentheses)."""
    rows, prev = [], None
    lev = sorted(summary["three_d"].items(), key=lambda kv: int(kv[0][2:]))
    for name, r in lev:
        cells = [name[2:], r.get("n_cells", "--")]
        for k in keys:
            x = r[group].get(k)
            o = order(prev[group].get(k), x, prev["h"], r["h"]) if prev else None
            cells.append(e(x) + (f" ({0.0 if abs(o) < 0.05 else o:.1f})" if o is not None else ""))
        for k, f in extra:
            cells.append(f(r))
        rows.append(cells)
        prev = r
    return rows


def wall(r):
    w = r.get("wall_time")
    if not w:
        return "--"
    return f"{w:.0f} s" if w < 90 else (f"{w / 60:.0f} min" if w < 5400 else f"{w / 3600:.1f} h")


# --------------------------------------------------------------------------- #
#  Helical-fiber analyses from the run folders (mid-length section)
# --------------------------------------------------------------------------- #

def zero_d_ref(case):
    """0D Richardson reference (cached by run_study.run_0d) at the end of the ramp."""
    import cylinder_3d as C3
    import run_study as RS
    root = os.path.join(HERE, case)
    m0 = C3.read_0d(os.path.join(root, "chamber_cylinder_0d.json"))
    ne = int(m0["values"]["num_elements"])
    out = os.path.join(root, "output")
    ref = RS.richardson(RS.run_0d(m0["config"], m0["name"], ne, out),
                        RS.run_0d(m0["config"], m0["name"], 2 * ne, out))
    T, _ = C3.pressure_ramp(m0["t"], m0["P"])
    return ref, int(np.argmin(np.abs(ref["t"] - T)))


def signed_mid(case, run, cache):
    """Endo and epi rho (3D/0D - 1) and the twist-rate ratio at mid-length,
    end of the ramp."""
    import cylinder_3d as C3
    key = f"{case}/{run}"
    if key in cache:
        return cache[key]
    d = os.path.join(HERE, case, "runs", run)
    if not os.path.exists(os.path.join(d, "results")):
        return None
    with open(os.path.join(d, "case.json")) as f:
        n_steps = json.load(f)["n_steps"]
    if len([x for x in os.listdir(os.path.join(d, "results")) if x.endswith(".vtu")]) != n_steps:
        return None
    ref, i = zero_d_ref(case)
    s = C3.post(d, fields=False, mid_section=True)
    rho, R = np.array(s["rho"][-1]), np.array(s["R"])
    rr = np.interp(R, ref["R"], ref["rho"][i])
    cache[key] = dict(endo=float(rho[0] / rr[0] - 1), epi=float(rho[-1] / rr[-1] - 1),
                      beta=float(s["beta"][-1] / ref["beta"][i]))
    return cache[key]


HELICAL = [("passive_helical_L2", r"$\pm60^\circ$ (test case)"),
           ("passive_helical_a30_L2", r"uniform $30^\circ$"),
           ("passive_helical_a60to0_L2", r"$+60^\circ \to 0^\circ$")]


def helical_signed_table(cache):
    rows = []
    for case, helix in HELICAL:
        for nr in (1, 2, 4, 8):
            m = signed_mid(case, f"struct-p2_nr{nr}", cache)
            if m:
                rows.append([helix, nr, f"{m['endo']:+.2e}", f"{m['epi']:+.2e}", f"{m['beta']:.4f}"])
    return table(["helix", "nr", r"endo $\rho$: 3D/0D$-1$", r"epi $\rho$: 3D/0D$-1$", r"twist rate 3D/0D"],
                 rows or [[PENDING] + [""] * 4],
                 "Helical fibers, $2L_0$ tube: signed differences at the mid-length section at "
                 "120~mmHg.", "tab:helical-signed")


# --------------------------------------------------------------------------- #
#  Tables from the study summaries
# --------------------------------------------------------------------------- #

KIN = [("V", "$V$"), ("rho", r"$\rho(R)$"), ("h", "$h$"), ("eta", r"$\eta$")]
KIN_HELIX = KIN + [("beta", r"$\beta$"), ("phi", r"$\varphi$")]
KIN_CLOSED = KIN + [("eps", r"$\varepsilon$")]
FIELDS = [("E_RR", "$E_{RR}$"), ("E_TT", r"$E_{\Theta\Theta}$"), ("J", "$J$"),
          ("lambda_f", r"$\lambda_f$"), ("sigma_vm", r"$\sigma_\mathrm{vm}$"),
          ("sigma_dev", r"$\sigma_{\theta\theta}-\sigma_{rr}$")]
ACTIVE = [("V", "$V$"), ("P", "$P$"), ("rho", r"$\rho(R)$"), ("SV", "SV"), ("W", "SW"), ("EF", "EF"),
          ("Q_in", r"$Q_\mathrm{in}$"), ("Q_out", r"$Q_\mathrm{out}$"), ("dPdt+", r"$\dot P_\mathrm{max}$")]


def passive_tables(case, title, kin=KIN):
    s = load(f"{case}/output/summary.json")
    if not s:
        return f"\\paragraph{{{title}}} {PENDING}\n"
    where = " (mid-length section)" if s.get("section") == "mid-length" else ""
    ref = s.get("reference", {}).get("name", "")
    t1 = table(["nr", "cells"] + [x for _, x in kin] + ["run time"],
               level_rows(s, [k for k, _ in kin], extra=[("wall", wall)]),
               f"{title}: maximum relative differences of the kinematics over the pressure ramp, "
               f"3D vs 0D ({ref}){where}; observed orders in parentheses.", f"tab:{case}-kin")
    t2 = table(["nr", "cells"] + [x for _, x in FIELDS],
               level_rows(s, [k for k, _ in FIELDS], group="field_errors"),
               f"{title}: strains (relative to the largest radial or hoop strain), $J-1$, "
               f"$\\lambda_f-1$ and Cauchy stresses (relative to their maximum){where}.",
               f"tab:{case}-fields")
    return t1 + "\n\n" + t2


def active_tables(case, title):
    s = load(f"{case}/output/summary.json")
    if not s:
        return f"\\paragraph{{{title}}} {PENDING}\n"
    keys = ACTIVE + ([("beta", r"$\beta$")] if s.get("twisting") else [])
    rows = level_rows(s, [k for k, _ in keys])
    for row, (name, r) in zip(rows, sorted(s["three_d"].items(), key=lambda kv: int(kv[0][2:]))):
        row.append(f"{r['errors'].get('valves', float('nan')):.1f}")
    for row in rows:
        row.pop(1)  # no cell counts in the active summaries
    return table(["nr"] + [x for _, x in keys] + ["valves [ms]"], rows,
                 f"{title}: maximum relative differences over cycle {s.get('cycle', 1)} "
                 "(flows relative to their peak; SV stroke volume, SW stroke work, EF ejection "
                 "fraction; valves: largest offset of the valve opening and closing times).",
                 f"tab:{case}")


def zero_d_table(cases):
    """cases: (case folder, the 0D solution used as the 3D reference, label or None)."""
    rows, seen = [], set()
    for case, used, name_ in cases:
        z = load(f"{case}/output/zero_d_convergence.json")
        if not z:
            continue
        best = z.get("reference", "").replace("richardson", "")
        for k, d in z["differences"].items():
            name = (f"Richardson {k[10:]}" if k.startswith("richardson") else f"ne = {k[2:]}")
            if k == used:
                name = rf"\textbf{{{name}}} (3D reference)"
            label = (name_ or case.replace("_", r"\_")) + f" ({best})"
            rows.append([label if label not in seen else "", name,
                         e(d.get("V")), e(d.get("P")), e(d.get("rho")), e(d.get("beta"))])
            seen.add(label)
    return table(["case (best estimate)", "0D solution", "$V$", "$P$", r"$\rho(R)$", r"$\beta$"],
                 rows or [[PENDING, "", "", "", "", ""]],
                 "0D block: maximum relative differences of direct solutions with ne radial "
                 "elements and of Richardson extrapolations (ne, 2ne) to the best estimate of the "
                 "converged 0D solution, the Richardson extrapolation of the two finest levels "
                 "(in parentheses). Bold: the reference of the 3D comparisons.", "tab:zero-d")


def direct_reference_table(cases):
    """3D errors against the Richardson reference and against the direct fine
    0D solution (summary_ne<N>.json of --reference-ne)."""
    import glob
    rows = []
    for case, keys in cases:
        s = load(f"{case}/output/summary.json")
        for path in sorted(glob.glob(os.path.join(HERE, case, "output", "summary_ne*.json"))):
            d = load(os.path.relpath(path, HERE))
            ne = os.path.basename(path)[len("summary_ne"):-len(".json")]
            if not s or not d:
                continue
            for lev in sorted(d["three_d"], key=lambda k: int(k[2:])):
                if lev not in s["three_d"]:
                    continue
                a, b = s["three_d"][lev]["errors"], d["three_d"][lev]["errors"]
                rows.append([case.replace("_", r"\_"), lev[2:], ne]
                            + [f"{e(a.get(k))} / {e(b.get(k))}" for k in keys])
    return table(["case", "nr", "ne", "$V$", r"$\rho(R)$", r"$\beta$ or $P$"],
                 rows or [[PENDING, "", "", "", "", ""]],
                 "3D vs 0D with the Richardson reference / with the direct 0D solution with ne "
                 "radial elements (no extrapolation): maximum relative differences.",
                 "tab:direct-ref")


def replay_table(cases):
    """0D loaded with the 3D cavity pressure (zero_d_replay.py)."""
    rows = []
    for case in cases:
        r = load(f"{case}/output/replay.json")
        if not r:
            continue
        for lev in sorted(r["levels"], key=lambda k: int(k[2:])):
            x = r["levels"][lev]
            rows.append([case.replace("_", r"\_"), lev[2:]]
                        + [e(x.get(k)) for k in ("V", "SV", "V_mid", "SV_mid", "rho", "beta")])
    return table(["case", "nr", "$V$ (tube)", "SV (tube)", "$V$ (mid)", "SV (mid)", r"$\rho(R)$ (mid)",
                  r"$\beta$ (mid)"], rows or [[PENDING] + [""] * 7],
                 "0D block loaded with the cavity pressure of the coupled 3D run: maximum relative "
                 "differences over the first cycle. Tube: the volume of the whole 3D tube; mid: "
                 "$V = \\pi(R_i + \\rho(R_i))^2 L$ and the wall state of the mid-length section "
                 "(eighth model: all nodes).", "tab:replay")


CLOSED_LENGTHS = [("passive_isotropic_closed_L2", 2), ("passive_isotropic_closed_L4", 4),
                  ("passive_isotropic_closed", 6), ("passive_isotropic_closed_L8", 8)]


def closed_length_table():
    """Free-length isotropic tube: differences at the mid-length plane against
    the tube length."""
    rows = []
    for case, n in CLOSED_LENGTHS:
        s = load(f"{case}/output/summary.json")
        if not s:
            continue
        for lev in sorted(s["three_d"], key=lambda k: int(k[2:])):
            r = s["three_d"][lev]
            rows.append([f"${n}L_0$" + (" (used)" if n == 6 else ""), lev[2:], r.get("n_cells", "--")]
                        + [e(r["errors"].get(k)) for k in ("V", "rho", "h", "eps")] + [wall(r)])
    return table(["tube", "nr", "cells", "$V$", r"$\rho(R)$", "$h$", r"$\varepsilon$", "run time"],
                 rows or [[PENDING] + [""] * 7],
                 "Free length, isotropic wall: maximum relative differences at the mid-length plane "
                 "over the pressure ramp against the length of the tube.", "tab:closed-length")


FIGURES = [
    ("passive_anisotropic_full/output/comparison_struct-p2.pdf",
     "Full Holzapfel-Ogden wall (eighth model): 0D vs 3D through the wall and mesh convergence."),
    ("passive_isotropic_closed/output/closed_end.pdf",
     "Free length: the 3D tube closed by an end plate. Near the plate the wall is held in; "
     "(a, b) the end layer against the tube length, (c) against the mesh, (d) the axial strain "
     "at the mid-length plane. The axial strain of the isotropic wall is small (0.5\\% at "
     "120~mmHg), so its relative differences in (b) are large.\\label{fig:closed-end}"),
    ("passive_helical_L2/output/helical_bcs_sketch.pdf",
     "Helical fibers: the uniform 0D state carries the section shear $\\sigma_{\\theta z}(R)$, "
     "which the free 3D ends cannot."),
    ("passive_helical_L2/output/helical_discrepancy.pdf",
     "Helical fibers ($\\pm60^\\circ$, $2L_0$ tube): radial expansion and twist along the axis, "
     "0D vs 3D. The end layers deviate; the middle follows the 0D block."),
    ("passive_helical_L2/output/comparison_struct-p2.pdf",
     "Helical fibers ($\\pm60^\\circ$, $2L_0$ tube): 0D vs 3D at the mid-length section."),
    ("active_prescribed/output/comparison_cycle1.pdf",
     "Active contraction, circumferential fibers: 0D vs coupled 3D, first cycle."),
    ("active_helical_L2/output/comparison_cycle1.pdf",
     "Active contraction, helical fibers ($2L_0$ tube): 0D vs coupled 3D, first cycle."),
]


def overview_figure(path):
    """One page with the convergence of all studies (log-log, mesh size h)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(15.0, 9.4))
    style = dict(lw=1.5, ms=5, mfc="none")
    H0 = 13.93  # wall thickness of the base geometry [mm]

    def levels(case):
        s = load(f"{case}/output/summary.json")
        if not s:
            return []
        return [dict(r, nr=int(k[2:])) for k, r in
                sorted(s["three_d"].items(), key=lambda kv: -kv[1]["h"])]

    def curve(a, case, key, label, c, m, ls="-", by_nr=False):
        """by_nr: plot at the mesh size of the base wall with the same nr
        (the wall-thickness variants)."""
        lev = [r for r in levels(case) if key in r["errors"]]
        if lev:
            x = [H0 / r["nr"] if by_nr else r["h"] * 1e3 for r in lev]
            a.loglog(x, [r["errors"][key] for r in lev], ls, marker=m,
                     color=c, label=label, **style)

    def slope(a, h0, e0, p):
        h = np.array([h0, h0 / 4])
        a.loglog(h, e0 * (h / h0) ** p, ":", color="0.5", lw=1)
        a.text(h[1], e0 * (h[1] / h0) ** p, f" $h^{p}$", color="0.4", va="top", fontsize=8)

    walls = [("passive_isotropic", "isotropic"), ("passive_anisotropic", "one fiber family"),
             ("passive_anisotropic_full", "full Holzapfel-Ogden")]
    for a, key, name, p in ((axs[0, 0], "V", "volume $V$", 4), (axs[0, 1], "rho", r"radial displacement $\rho(R)$", 3)):
        for i, (case, label) in enumerate(walls):
            curve(a, case, key, label, f"C{i}", "osd"[i])
        for case, label, m in (("passive_full_thin", r"full HO, $R_e/R_i = 1.2$ (at its nr)", "<"),
                               ("passive_full_thick", r"full HO, $R_e/R_i = 2.5$ (at its nr)", ">")):
            curve(a, case, key, label, "C2", m, "--", by_nr=True)
        curve(a, "passive_isotropic_soft", key, r"isotropic, $\mu/2$ (large deformation)", "C0", "x", "--")
        for i, (case, _) in enumerate(walls):
            curve(a, case + "_closed", key, None, f"C{i}", "osd"[i], "-.")
        a.plot([], [], "k-.", label="free length (closed ends, mid-length plane)")
        slope(a, 14.0, 2e-3 if key == "V" else 1e-2, p)
        a.set(title=f"(a{1 if key == 'V' else 2}) passive, circumferential fibers: {name}")
    for a, key, name in ((axs[0, 2], "rho", r"$\rho(R)$"), (axs[1, 0], "beta", r"twist rate $\beta$")):
        for i, (case, label) in enumerate(HELICAL):
            curve(a, case, key, label, f"C{i + 3}", "o^s"[i])
        a.set(title=f"({'b1' if key == 'rho' else 'b2'}) passive, helical fibers, $2L_0$ tube, "
                    f"mid-length: {name}")
    a = axs[1, 1]
    for i, case in enumerate(["passive_anisotropic", "passive_anisotropic_full", "passive_helical",
                              "active_prescribed"]):
        z = load(f"{case}/output/zero_d_convergence.json")
        if not z:
            continue
        d = z["differences"]
        ne = sorted(int(k[2:]) for k in d if k.startswith("ne"))
        a.loglog(ne, [d[f"ne{n}"]["V"] for n in ne], "-", marker="osd^"[i], color=f"C{i}",
                 label=case.replace("passive_", "").replace("_", " ").replace("helical", r"helical $\pm60^\circ$"),
                 **style)
        rich = sorted(int(k[10:].split("/")[0]) for k in d if k.startswith("richardson"))
        a.loglog(rich, [d[f"richardson{n}/{2 * n}"]["V"] for n in rich], "--",
                 marker="osd^"[i], color=f"C{i}", **style)
    a.plot([], [], "k-", label="direct, ne elements")
    a.plot([], [], "k--", label="Richardson ne/2ne")
    a.set(title="(c) 0D block: volume vs its best estimate", xlabel="0D radial elements ne (Richardson: ne/2ne)")
    a.set_xticks([16, 32, 64, 128, 256], ["16", "32", "64", "128", "256"])
    a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    a = axs[1, 2]
    active = [("active_prescribed", "circ., 400 steps"), ("active_prescribed_dt2", "circ., 800 steps"),
              ("active_helical_L2", r"helical, $2L_0$ tube")]
    for i, (case, label) in enumerate(active):
        curve(a, case, "P", f"$P$, {label}", f"C{i}", "osd^"[i])
        curve(a, case, "V", f"$V$, {label}", f"C{i}", "osd^"[i], "--")
    curve(a, "active_thick", "P", r"$P$, circ., $R_e/R_i = 2.5$ (at its nr)", "C4", "v", by_nr=True)
    curve(a, "active_thick", "V", r"$V$, circ., $R_e/R_i = 2.5$ (at its nr)", "C4", "v", "--", by_nr=True)
    a.set(title="(d) active, coupled: pressure (solid), volume (dashed)")
    for a in axs.flat:
        if a is not axs[1, 1]:
            a.set(xlabel="mesh size $h = (R_e - R_i)$/nr [mm]")
            a.set_xticks([0.875, 1.75, 3.5, 7.0, 14.0], ["0.875\n(nr 16)", "1.75 (8)", "3.5 (4)", "7 (2)", "14 (1)"])
            a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            a.set_xlim(0.75, 17.0)
        a.set(ylabel="max relative difference")
        a.grid(True, which="major", color="0.9", lw=0.6)
        a.legend(fontsize=7.5)
        a.title.set_fontsize(9.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def figures():
    out = []
    for path, cap in FIGURES:
        full = os.path.join(HERE, path)
        if os.path.exists(full):
            out.append(r"\begin{figure}[p]\centering"
                       rf"\includegraphics[width=\textwidth,height=0.85\textheight,keepaspectratio]{{{full}}}"
                       rf"\caption{{{cap}}}\end{{figure}}")
    return "\n".join(out)


# --------------------------------------------------------------------------- #
#  Document
# --------------------------------------------------------------------------- #

def read_text(name):
    """Hand-written narrative sections (summary/text/<name>.tex), if present."""
    path = os.path.join(OUT, "text", f"{name}.tex")
    if os.path.exists(path):
        with open(path) as f:
            return f.read()
    return ""


PREAMBLE = r"""\documentclass[10pt]{article}
\usepackage[margin=2cm]{geometry}
\usepackage{booktabs,graphicx,amsmath,amssymb,float}
\usepackage[hidelinks]{hyperref}
\usepackage[font=small,labelfont=bf]{caption}
\setlength{\parindent}{0pt}\setlength{\parskip}{4pt}
\title{ChamberCylinder: verification of the 0D block against 3D finite elements}
\author{svZeroDSolver, branch \texttt{claude/chamber-cylinder-3d-verification}}
\date{\today}
\begin{document}
\maketitle
\tableofcontents
"""


def build():
    cache_path = os.path.join(OUT, "helical_cache.json")
    cache = load(os.path.relpath(cache_path, HERE)) or {}
    parts = [PREAMBLE]
    overview = overview_figure(os.path.join(OUT, "convergence_overview.pdf"))
    parts.append(r"\section{Summary}" + "\n" + read_text("summary")
                 + "\n" + r"\begin{figure}[p]\centering"
                 + rf"\includegraphics[width=\textwidth]{{{overview}}}"
                 + r"\caption{Convergence of all studies: maximum relative differences between the 3D "
                 + r"and the 0D solutions against the mesh size (a, b, d) and of the 0D block against "
                 + r"its converged solution (c). Helical fibers: the mid-length section.}"
                 + r"\label{fig:overview}\end{figure}")
    parts.append(r"\section{Setup}" + "\n" + read_text("setup"))
    parts.append(r"\section{Passive inflation, circumferential fibers (eighth model)}" + "\n"
                 + read_text("passive"))
    for case, title in (("passive_isotropic", "Isotropic neo-Hookean wall"),
                        ("passive_anisotropic", "One circumferential fiber family"),
                        ("passive_anisotropic_full", "Full Holzapfel-Ogden wall")):
        parts.append(rf"\subsection{{{title} (\texttt{{{case.replace('_', chr(92) + '_')}}})}}" + "\n"
                     + passive_tables(case, title))
    parts.append(r"\subsection{Wall thickness}\label{sec:wall}" + "\n" + read_text("wall")
                 + passive_tables("passive_full_thin", r"Thin wall, $R_e/R_i = 1.2$")
                 + "\n\n" + passive_tables("passive_full_thick", r"Very thick wall, $R_e/R_i = 2.5$"))
    parts.append(r"\subsection{Large deformation}\label{sec:large}" + "\n" + read_text("large")
                 + passive_tables("passive_isotropic_soft", r"Isotropic wall, $\mu = 50$~kPa"))
    parts.append(r"\subsection{Free length (closed ends)}\label{sec:closed}" + "\n" + read_text("closed"))
    for case, title in (("passive_isotropic_closed", "Isotropic wall, free length"),
                        ("passive_anisotropic_closed", "One circumferential fiber family, free length"),
                        ("passive_anisotropic_full_closed", "Full Holzapfel-Ogden wall, free length")):
        parts.append(rf"\subsubsection{{{title}}}" + "\n" + passive_tables(case, title, KIN_CLOSED))
    parts.append(r"\subsubsection{Tube length}\label{sec:closed-length}" + "\n" + read_text("closed_length") + "\n"
                 + closed_length_table())
    parts.append(r"\section{Passive inflation, helical fibers (full cylinder, $2L_0$)}\label{sec:helical}"
                 + "\n" + read_text("helical"))
    parts.append(r"\subsection{Helix of the test case, $\pm60^\circ$}" + "\n"
                 + passive_tables("passive_helical_L2", r"Helical fibers, $\pm60^\circ$, $2L_0$", KIN_HELIX))
    parts.append(r"\subsection{Uniform $30^\circ$ helix (large twist)}" + "\n"
                 + read_text("helical_a30")
                 + passive_tables("passive_helical_a30_L2", r"Uniform $30^\circ$ helix, $2L_0$", KIN_HELIX))
    parts.append(r"\subsection{Helix from $+60^\circ$ to $0^\circ$ (varying angle, large twist)}"
                 + "\n" + read_text("helical_a60to0")
                 + passive_tables("passive_helical_a60to0_L2", r"Helix $+60^\circ$ to $0^\circ$, $2L_0$",
                                  KIN_HELIX))
    parts.append(r"\subsection{Signed differences}" + "\n" + helical_signed_table(cache))
    parts.append(r"\section{Accuracy of the 0D reference}\label{sec:zero-d}" + "\n" + read_text("zero_d")
                 + zero_d_table([("passive_anisotropic", "richardson32/64", None),
                                 ("passive_anisotropic_full", "richardson32/64", None),
                                 ("passive_full_thin", "richardson32/64", None),
                                 ("passive_full_thick", "richardson32/64", None),
                                 ("passive_isotropic_soft", "richardson32/64", None),
                                 ("passive_helical", "richardson32/64", r"helical $\pm60^\circ$"),
                                 ("active_prescribed", "richardson16/32", None)])
                 + "\n\n" + read_text("zero_d_direct")
                 + direct_reference_table([("passive_anisotropic", ("V", "rho", "beta")),
                                           ("passive_anisotropic_full", ("V", "rho", "beta")),
                                           ("active_prescribed", ("V", "rho", "P"))]))
    parts.append(r"\section{Active contraction coupled to a circulation}" + "\n" + read_text("active"))
    for case, title in (("active_prescribed", "Circumferential fibers, 400 steps per cycle"),
                        ("active_prescribed_dt2", "Circumferential fibers, 800 steps per cycle"),
                        ("active_prescribed_c2", "Circumferential fibers, second (periodic) cycle"),
                        ("active_thick", r"Circumferential fibers, thick wall ($R_e/R_i = 2.5$)"),
                        ("active_thick_dt2", r"Thick wall, 800 steps per cycle"),
                        ("active_helical_L2", r"Helical fibers ($\pm60^\circ$, $2L_0$ tube)")):
        parts.append(rf"\subsection{{{title}}}" + "\n" + active_tables(case, title))
    parts.append(r"\subsection{Diagnostic: the 0D block at the 3D pressure}\label{sec:replay}" + "\n"
                 + read_text("replay")
                 + replay_table(["active_prescribed", "active_helical_L2"]))
    parts.append(r"\section{Numerical notes}" + "\n" + read_text("notes"))
    parts.append(r"\section{Open points}" + "\n" + read_text("open"))
    parts.append(r"\clearpage\appendix\section{Figures}" + "\n" + figures())
    parts.append(r"\end{document}")
    os.makedirs(OUT, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(cache, f, indent=1)
    tex = os.path.join(OUT, "verification_summary.tex")
    with open(tex, "w") as f:
        f.write("\n\n".join(parts))
    return tex


def compile_pdf(tex):
    exe = shutil.which("pdflatex") or "/Library/TeX/texbin/pdflatex"
    for _ in range(2):  # table of contents
        r = subprocess.run([exe, "-interaction=nonstopmode", "-halt-on-error", os.path.basename(tex)],
                           cwd=os.path.dirname(tex), capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-3000:])
            raise SystemExit("pdflatex failed")
    return tex[:-4] + ".pdf"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-compile", action="store_true")
    a = ap.parse_args()
    tex = build()
    print(compile_pdf(tex) if not a.no_compile else tex)


if __name__ == "__main__":
    main()
