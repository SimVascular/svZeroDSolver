# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Passive inflation of a fixed-length cylinder: 0D vs. 3D mesh convergence.

CASE is a folder with the 0D input chamber_cylinder_0d.json, which defines the
problem. For each 3D mesh level (--nr, elements through the wall) the script
writes the svMultiPhysics case into CASE/runs/, runs it if --run is given
(Docker image simvascular/solver by default, or --exe), and compares over the
pressure ramp the kinematic quantities of the 0D block:

  V       cavity volume                        |V - V_ref| / V_ref
  rho     radial displacement rho(R)           max |rho - rho_ref| / max |rho_ref|
  h       wall-thickness change rho(Re)-rho(Ri) relative to its maximum
  beta    twist per unit length                relative to max |beta_ref|
  phi     transmural rotation phi(R)-phi(Ri)   relative to the twist max |beta_ref| L
  eta     axial displacement eta(R)-eta(Ri)    relative to max |rho_ref|
  eps     axial strain (free length)           relative to max |eps_ref|

The profiles are compared at the radial node levels of the 3D mesh, which are
nodes of the 0D reference. The reference is the Richardson extrapolation of the
0D block at num_elements and 2 num_elements (second order), or the closed form
of the incompressible neo-Hookean tube when the wall is isotropic and
incompressible. Each error is the maximum over the 0D time points of the ramp;
the 3D model takes --substeps load steps per 0D step.

Strains and stresses are compared through the wall as well: level means of the
3D nodal fields (cylinder_3d.struct_fields) vs the 0D profiles of the same
quantities (zero_d_fields), interpolated to the 3D levels:

  E_RR, E_TT, E_ZZ  Green-Lagrange strain     relative to max |E_RR|, |E_TT|
  J                 Jacobian, element means   relative to max |J - 1|
  lambda_f          fiber stretch             relative to max |lambda_f - 1|
  sigma_vm          von Mises stress          relative to its maximum
  sigma_dev         sigma_tt - sigma_rr       relative to its maximum
  sigma_rr          radial Cauchy stress, element means, relative to P_max
  sigma_rr_endo     sigma_rr(Ri) + P (traction boundary condition), relative to P_max

The 3D hydrostatic stress oscillates between the node levels of an element,
as does J (nodal pressure oscillation of the displacement formulation), so J
and sigma_rr are compared as element means and the other stresses are
deviatoric.

Run from the repository root after building the Python module:
  PYTHONPATH=build/python python scripts/ChamberCylinder/verification_3d/run_study.py \\
      scripts/ChamberCylinder/verification_3d/passive_isotropic --nr 1 2 4 [--run] [--np 4]
Results: CASE/runs/<formulation>_nr<N>/ and CASE/output/ (summary.json,
comparison_<formulation>.pdf).
Pressures are reported in mmHg; the solver inputs are SI.
"""
import argparse
import copy
import glob
import hashlib
import json
import math
import os
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import cylinder_3d as C3  # noqa: E402
from ChamberCylinder_verification import fixed_length_exact  # noqa: E402

PROFILES = ("rho", "phi", "eta")
QUANTITIES = [  # key, label, plot marker
    ("V", "volume", "s-"), ("rho", "radial displacement", "o-"),
    ("h", "wall-thickness change", "D-"), ("beta", "twist", "^-"),
    ("phi", "transmural rotation", "v-"), ("eta", "axial displacement", "x-"),
    ("eps", "axial strain", "<-")]
FIELD_QUANTITIES = [  # key, label, plot marker
    ("E_RR", "radial strain", "o-"), ("E_TT", "hoop strain", "s-"), ("E_ZZ", "axial strain", "^-"),
    ("J", "Jacobian, element mean", "D-"), ("lambda_f", "fiber stretch", "v-"),
    ("sigma_vm", "von Mises stress", "x-"), ("sigma_dev", "sigma_tt - sigma_rr", "+-"),
    ("sigma_rr", "sigma_rr, element mean", "*-"), ("sigma_rr_endo", "sigma_rr(Ri) + P", "p-")]


def zero_d_fields(R, rho, phi, eta, beta, values, alpha=(0.0, 0.0), eps=0.0, stresses=True):
    """Through-wall profiles of strain, Jacobian, fiber stretch and Cauchy stress
    of a 0D solution.

    rho, phi, eta are nodal fields (n_t, n_R) at the radii R and beta the twist
    per unit length (n_t). With r = R + rho, the deformation gradient in the
    reference (e_R, e_Th, e_Z) and deformed (e_r, e_th, e_z) bases is
    [[1 + rho', 0, 0], [r phi', r/R, r beta], [eta', 0, 1 + eps]], with the
    radial derivatives from cubic splines of the nodal values. Strains are
    Green-Lagrange components in the reference basis, stresses Cauchy
    components in the deformed basis from the 0D material law (check_materials.
    stress_0d, the ST91 pressure of the local J). alpha: helix angles at R[0],
    R[-1] in degrees (fiber stretch). eps: axial strain, a number or one per
    time (n_t).
    """
    from scipy.interpolate import CubicSpline
    from check_materials import stress_0d
    R = np.asarray(R, float)
    rho, phi, eta = (np.atleast_2d(np.asarray(x, float)) for x in (rho, phi, eta))
    beta = np.atleast_1d(np.asarray(beta, float))
    eps = np.broadcast_to(np.asarray(eps, float), beta.shape)
    a = np.radians(alpha[0] + (alpha[1] - alpha[0]) * (R - R[0]) / (R[-1] - R[0]))
    vals = {k: float(values.get(k, 0.0)) for k in ("C1", "C2", "C3", "C4", "C5", "C6")}
    vals.update({k: float(x) for k, x in values.items()
                 if k.startswith("ho_") or k == "bulk_modulus"})
    keys = ["E_RR", "E_TT", "E_ZZ", "J", "lambda_f"]
    if stresses:
        keys += ["sigma_rr", "sigma_tt", "sigma_zz", "sigma_vm"]
    out = {k: np.zeros((len(beta), len(R))) for k in keys}
    s = np.array([1.0, 0.0, 0.0])
    for i in range(len(beta)):
        d = {k: CubicSpline(R, x[i])(R, 1) for k, x in (("rho", rho), ("phi", phi), ("eta", eta))}
        for j, Rj in enumerate(R):
            r = Rj + rho[i, j]
            F = np.array([[1.0 + d["rho"][j], 0.0, 0.0],
                          [r * d["phi"][j], r / Rj, r * beta[i]],
                          [d["eta"][j], 0.0, 1.0 + eps[i]]])
            C = F.T @ F
            f = np.array([0.0, np.cos(a[j]), np.sin(a[j])])
            J = np.linalg.det(F)
            for k, c in (("E_RR", 0), ("E_TT", 1), ("E_ZZ", 2)):
                out[k][i, j] = 0.5 * (C[c, c] - 1.0)
            out["J"][i, j] = J
            out["lambda_f"][i, j] = np.sqrt(f @ C @ f)
            if stresses:
                sig = F @ stress_0d(F, f, s, vals) @ F.T / J
                dev = sig - np.trace(sig) * np.eye(3) / 3.0
                for k, c in (("sigma_rr", 0), ("sigma_tt", 1), ("sigma_zz", 2)):
                    out[k][i, j] = sig[c, c]
                out["sigma_vm"][i, j] = np.sqrt(1.5 * np.sum(dev * dev))
    return out


def element_means(prof, order):
    """Mean of nodal profiles (n_t, n_R) over each element of a mesh with
    `order` + 1 equidistant node levels per element (Simpson for order 2,
    trapezoid for order 1)."""
    prof = np.asarray(prof)
    w = np.array([1.0, 4.0, 1.0]) / 6.0 if order == 2 else np.array([0.5, 0.5])
    n_el = (prof.shape[1] - 1) // order
    return np.stack([prof[:, order * e:order * e + order + 1] @ w for e in range(n_el)], axis=1)


def field_errors(sol, ref_fields, ref_R, P, order):
    """Max relative differences over the ramp of the field profiles of sol
    (sol["fields"] at the levels sol["R"], elements of `order` + 1 levels) vs
    the 0D reference profiles. J and sigma_rr, which carry the hydrostatic
    pressure, are compared as element means."""
    from scipy.interpolate import CubicSpline
    R = np.asarray(sol["R"])
    f3 = {k: np.asarray(v) for k, v in sol["fields"].items()}
    spline = lambda k: [CubicSpline(ref_R, row) for row in ref_fields[k]]  # noqa: E731
    at = lambda k: np.array([s(R) for s in spline(k)])  # noqa: E731
    rel = lambda k, scale: float(np.max(np.abs(f3[k] - at(k))) / scale)  # noqa: E731
    edges = R[::order]

    def rel_mean(k, scale):  # element means of the 3D levels and of the 0D spline
        ref = np.array([[s.integrate(a, b) / (b - a) for a, b in zip(edges[:-1], edges[1:])]
                        for s in spline(k)])
        return float(np.max(np.abs(element_means(f3[k], order) - ref)) / scale)
    sE = max(np.max(np.abs(ref_fields["E_RR"])), np.max(np.abs(ref_fields["E_TT"])))
    e = {k: rel(k, sE) for k in ("E_RR", "E_TT", "E_ZZ")}
    e["J"] = rel_mean("J", np.max(np.abs(ref_fields["J"] - 1.0)))
    e["lambda_f"] = rel("lambda_f", np.max(np.abs(ref_fields["lambda_f"] - 1.0)))
    if "sigma_rr" in ref_fields and "sigma_rr" in f3:
        Pm = np.max(np.abs(P))
        e["sigma_vm"] = rel("sigma_vm", np.max(ref_fields["sigma_vm"]))
        dev_ref = at("sigma_tt") - at("sigma_rr")
        e["sigma_dev"] = float(np.max(np.abs(f3["sigma_tt"] - f3["sigma_rr"] - dev_ref))
                               / np.max(np.abs(dev_ref)))
        e["sigma_rr"] = rel_mean("sigma_rr", Pm)
        e["sigma_rr_endo"] = float(np.max(np.abs(f3["sigma_rr"][:, 0] + np.asarray(P))) / Pm)
    return e


def run_0d(config, name, num_elements=None, cache_dir=None):
    """0D pressure, volume, twist and nodal profiles rho, phi, eta over time.

    Results are cached in cache_dir under a hash of the 0D input.
    """
    cfg = copy.deepcopy(config)
    blk = next(v for v in cfg["vessels"] if v["vessel_name"] == name)
    vals = blk["zero_d_element_values"]
    if num_elements:
        vals["num_elements"] = num_elements
    ne = int(vals.get("num_elements", 10))
    cache = None
    if cache_dir:
        key = hashlib.sha1(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16]
        cache = os.path.join(cache_dir, f"zero_d_{ne}_{key}.npz")
        if os.path.exists(cache):
            with np.load(cache) as d:
                return {k: d[k] for k in d.files}
    import svzerod
    tic = time.time()
    out = svzerod.simulate(cfg)

    def series(var):
        d = out[out["name"] == var].sort_values("time")
        return d["time"].to_numpy(), d["y"].to_numpy()
    t, V = series(f"volume:{name}")
    pname = next(n for n in out["name"].unique() if n.startswith("pressure:") and name in n)
    sol = dict(t=t, P=series(pname)[1], V=V, beta=series(f"beta:{name}")[1],
               eps=series(f"eps:{name}")[1], R=np.linspace(vals["Ri"], vals["Re"], ne + 1))
    for f in PROFILES:
        sol[f] = np.stack([series(f"{f}_{A}:{name}")[1] for A in range(ne + 1)], axis=1)
    if cache:
        os.makedirs(cache_dir, exist_ok=True)
        np.savez(cache, **sol)
    print(f"0D ne = {ne}: {time.time() - tic:.0f} s", flush=True)
    return sol


def exact_isotropic(t, P, Ri, Re, L, mu, n=257):
    """Closed form of the fixed-length incompressible neo-Hookean tube."""
    R = np.linspace(Ri, Re, n)
    ri = np.array([fixed_length_exact(p, Ri, Re, mu) for p in P])
    rho = np.sqrt(R[None, :] ** 2 + ri[:, None] ** 2 - Ri ** 2) - R[None, :]
    z = np.zeros_like(rho)
    return dict(t=t, P=P, V=math.pi * ri ** 2 * L, beta=np.zeros_like(P), eps=np.zeros_like(P),
                R=R, rho=rho, phi=z, eta=z)


def richardson(coarse, fine):
    """Second-order Richardson extrapolation of two 0D solutions (ne, 2 ne) on
    the coarse nodes."""
    ext = dict(coarse)
    for k in ("V", "beta", "eps"):
        if k in coarse and k in fine:  # older cached solutions have no eps
            ext[k] = (4.0 * fine[k] - coarse[k]) / 3.0
    for k in PROFILES:
        ext[k] = (4.0 * fine[k][:, ::2] - coarse[k]) / 3.0
    return ext


def restrict(sol, t):
    """Solution at the times t (a subset of its time points)."""
    idx = [int(np.argmin(np.abs(sol["t"] - ti))) for ti in t]
    if not np.allclose(sol["t"][idx], t, rtol=0, atol=1e-9 * max(t)):
        raise ValueError("the compared solutions have different time points")
    pick = lambda v: {k: pick(w) for k, w in v.items()} if isinstance(v, dict) \
        else np.asarray(v)[idx]  # noqa: E731
    return {k: (v if k == "R" or np.ndim(v) == 0 else pick(v)) for k, v in sol.items()}


def tolist(x):
    """JSON-ready copy of an array or a dict of arrays."""
    return {k: tolist(v) for k, v in x.items()} if isinstance(x, dict) else np.asarray(x).tolist()


def has_twist(ref, L):
    """True when the reference twists beyond round-off (helical fibers). The
    0D Newton tolerance leaves end-face rotations of up to ~4e-7 rho without
    fibers (large deformation); the smallest helical twist is ~2e-2 rho."""
    twist = np.max(np.abs(ref["beta"])) * L * np.max(ref["R"])  # end-face rotation
    return bool(twist > 1e-5 * np.max(np.abs(ref["rho"])))


def errors(sol, ref, L):
    """Max relative differences over the ramp of sol vs ref (same times)."""
    at = lambda f: np.array([np.interp(sol["R"], ref["R"], row) for row in ref[f]])  # noqa: E731
    rho, rho_ref = np.asarray(sol["rho"]), at("rho")
    dh, dh_ref = rho[:, -1] - rho[:, 0], rho_ref[:, -1] - rho_ref[:, 0]
    scale_rho = np.max(np.abs(rho_ref))
    e = dict(V=np.max(np.abs(np.asarray(sol["V"]) / ref["V"] - 1.0)),
             rho=np.max(np.abs(rho - rho_ref)) / scale_rho,
             h=np.max(np.abs(dh - dh_ref)) / np.max(np.abs(dh_ref)),
             eta=np.max(np.abs(np.asarray(sol["eta"]) - at("eta"))) / scale_rho)
    e_ref = np.max(np.abs(ref.get("eps", 0.0)))
    if e_ref > 0.0:  # free length (closed ends)
        e["eps"] = np.max(np.abs(np.asarray(sol["eps"]) - ref["eps"])) / e_ref
    b_ref = np.max(np.abs(ref["beta"]))
    if has_twist(ref, L):  # helical fibers
        e["beta"] = np.max(np.abs(np.asarray(sol["beta"]) - ref["beta"])) / b_ref
        e["phi"] = np.max(np.abs(np.asarray(sol["phi"]) - at("phi"))) / (b_ref * L)
    return {k: float(v) for k, v in e.items()}


def table_row(label, e, quantities=QUANTITIES, w=12):
    return f"{label:<12s}" + "".join(f"{e[k]:>{w}.2e}" if k in e else f"{'-':>{w}s}"
                                     for k, _, _ in quantities)


def table_header(quantities=QUANTITIES, w=12):
    return f"{'':<12s}" + "".join(f"{k:>{w}s}" for k, _, _ in quantities)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case", help="folder with chamber_cylinder_0d.json")
    ap.add_argument("--nr", type=int, nargs="+", default=[1, 2, 4])
    ap.add_argument("--formulation", default="struct-p2", choices=sorted(C3.FORMULATIONS))
    ap.add_argument("--sector", type=float, help="default: 90, or 360 with helical fibers")
    ap.add_argument("--full-length", action="store_true",
                    help="mesh the full length (default: half, mid-length symmetry plane)")
    ap.add_argument("--run", action="store_true",
                    help="run svMultiPhysics where complete results of the same case are missing")
    ap.add_argument("--rerun", action="store_true", help="with --run: also rerun complete cases")
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--exe", help="local svmultiphysics executable instead of Docker")
    ap.add_argument("--max-time", type=float, default=1800.0,
                    help="wall-time limit per 3D run [s]")
    ap.add_argument("--substeps", type=int, default=1, help="3D load steps per 0D time step")
    ap.add_argument("--ls-type", default="BICG")
    ap.add_argument("--algebra", default="fsils")
    ap.add_argument("--preconditioner", default="fsils")
    ap.add_argument("--ls", choices=["trilinos-ml"],
                    help="preset linear solver (default: --ls-type with fsils)")
    ap.add_argument("--end-damping", type=float,
                    help="full cylinder: damping of the end dashpots in units of mu / Re * dt "
                         "other than the default 1e-3 (runs/<formulation>_nr<n>_damp<value>)")
    ap.add_argument("--reference-ne", type=int,
                    help="compare with the direct 0D solution with this many radial elements "
                         "instead of the Richardson extrapolation (output/summary_ne<N>.json, "
                         "no figure)")
    a = ap.parse_args()
    sys.stdout.reconfigure(line_buffering=True)

    case_root = os.path.abspath(a.case)
    model_0d = os.path.join(case_root, "chamber_cylinder_0d.json")
    m0 = C3.read_0d(model_0d)
    v = m0["values"]
    mat = C3.material(v, m0["formulation"])
    fibers = mat["model"] in ("HGO", "HolzapfelOgden")
    # free length (0D ends "lids", pressure-loaded lids): the 3D tube is
    # closed by an end plate (cylinder_3d.generate closed_end)
    closed = m0["formulation"]["ends"] == "lids"
    # helical fibers and closed ends: the 3D ends differ from the uniform 0D
    # state (end layers), so only the mid-length section is compared
    mid = C3.helical(mat) or closed
    T, Pmax = C3.pressure_ramp(m0["t"], m0["P"])
    Ri, Re, L = v["Ri"], v["Re"], v["length"]
    mmhg = C3.MMHG

    # time points of the 0D input on the ramp; the 3D case uses the same steps
    sp = m0["config"]["simulation_parameters"]
    t_all = np.linspace(0.0, sp["cardiac_period"], int(sp["number_of_time_pts_per_cardiac_cycle"]))
    t0 = t_all[t_all <= T * (1 + 1e-12)]
    out_dir = os.path.join(case_root, "output")
    zero_d = lambda ne=None: restrict(run_0d(m0["config"], m0["name"], ne, out_dir), t0)  # noqa: E731

    # 3D: generate and run each mesh level (stop at the first failure)
    ls = dict(type=a.ls_type, algebra=a.algebra, preconditioner=a.preconditioner)
    if a.ls == "trilinos-ml":  # GMRES with algebraic multigrid: needed from nr = 16 on
        ls = dict(type="GMRES", algebra="trilinos", preconditioner="trilinos-ml",
                  assembly="fsils", tol=1e-6, max_iter=1000, krylov=200)
    support, suffix = {}, ""
    if a.end_damping is not None:
        support, suffix = dict(end_damping=a.end_damping), f"_damp{a.end_damping:g}"
    cases = {}
    for nr in a.nr:
        case_dir = os.path.join(case_root, "runs", f"{a.formulation}_nr{nr}{suffix}")
        os.makedirs(case_dir, exist_ok=True)
        old = None
        if os.path.exists(os.path.join(case_dir, "case.json")):
            with open(os.path.join(case_dir, "case.json")) as f:
                old = json.load(f)
        if not a.run and old is not None:  # post-process the existing case as it was run
            case = old
        else:
            case = C3.generate(model_0d, case_dir, nr, a.formulation, a.sector,
                               n_steps=a.substeps * (len(t0) - 1), linear_solver=ls,
                               half_length=not a.full_length, closed_end=closed, **support)
        cases[nr] = (case_dir, case)
        done = len(glob.glob(os.path.join(case_dir, "results", "result_*.vtu"))) == case["n_steps"]
        if a.run and not a.rerun and done and old == json.loads(json.dumps(case)):
            print(f"3D nr = {nr}: complete results of the same case exist, not rerun")
        elif a.run:
            print(f"3D nr = {nr}: running ({case['n_cells']} tets, {case['n_steps']} steps)")
            shutil.rmtree(os.path.join(case_dir, "results"), ignore_errors=True)  # no stale steps
            try:
                wall = C3.run(case_dir, a.np, a.exe, max_time=a.max_time)
            except C3.SolverFailure as err:
                print(f"3D nr = {nr}: FAILED - {err}")
                break  # finer meshes would fail too
            with open(os.path.join(case_dir, "wall_time.txt"), "w") as f:
                f.write(f"{wall:.1f}\n")
            print(f"3D nr = {nr}: done in {wall:.0f} s")

    # 0D reference
    if a.reference_ne:
        ref, ref_name = zero_d(a.reference_ne), f"0D ne = {a.reference_ne}"
    elif mat["kappa"] > 0.0 or fibers:
        ne = int(v["num_elements"])
        ref, ref_name = richardson(zero_d(ne), zero_d(2 * ne)), f"0D Richardson ne = {ne}/{2 * ne}"
    else:
        ref, ref_name = exact_isotropic(t0, Pmax * t0 / T, Ri, Re, L, mat["mu"]), "exact"
    # strain and stress profiles of the 0D (exact: no stress, the pressure is not
    # a function of J)
    alpha = (mat.get("alpha_endo", 0.0), mat.get("alpha_epi", 0.0))
    fields_of = lambda s, stresses=True: zero_d_fields(  # noqa: E731
        s["R"], s["rho"], s["phi"], s["eta"], s["beta"], v, alpha, eps=s.get("eps", 0.0),
        stresses=stresses)
    ref_fields = fields_of(ref, ref_name != "exact")
    summary = dict(case=os.path.basename(case_root), material=mat, P_max_mmHg=Pmax / mmhg,
                   reference=dict(name=ref_name, fields=tolist(ref_fields),
                                  **{k: np.asarray(x).tolist() for k, x in ref.items()}),
                   section="mid-length" if mid else "all", zero_d={}, three_d={})
    twist = np.degrees(ref["beta"][-1] * L)
    print(f"0D reference ({ref_name}) at {Pmax / mmhg:.0f} mmHg: V = {ref['V'][-1] * 1e6:.4f} mL, "
          f"rho endo/epi = {ref['rho'][-1][0] * 1e3:.4f}/{ref['rho'][-1][-1] * 1e3:.4f} mm, "
          f"twist over the length = {twist:.4f} deg")
    if closed:
        print(f"free length: axial strain at {Pmax / mmhg:.0f} mmHg {ref['eps'][-1]:.5f}; 3D closed by "
              "an end plate, compared at the mid-length plane (V = pi (Ri + rho(Ri))^2 (1 + eps) L)")
    elif mid:
        print("helical fibers: 3D compared at the mid-length section (middle element layer; "
              "beta = local twist rate, V = pi (Ri + rho(Ri))^2 L)")

    # compare
    rows, field_rows = [], []
    for nr, (case_dir, case) in cases.items():
        try:
            out = C3.post(case_dir, mid_section=mid)
        except FileNotFoundError:
            print(f"3D nr = {nr}: generated {case_dir} ({case['n_cells']} tets), no results yet")
            continue
        sol = {k: (x if isinstance(x, dict) else np.asarray(x)) for k, x in out.items()}
        if len(sol["t"]) != case["n_steps"] + 1:
            print(f"3D nr = {nr}: incomplete results ({len(sol['t']) - 1} of {case['n_steps']} steps)")
            continue
        sol = restrict(sol, t0)
        e = errors(sol, ref, L)
        fe = field_errors(sol, ref_fields, ref["R"], ref["P"], case["order"])
        rec = dict(h=(Re - Ri) / nr, n_cells=case["n_cells"], n_nodes=case["n_nodes"], errors=e,
                   field_errors=fe, order=case["order"],
                   non_homogeneity=float(np.max(sol["rho_spread"]) / np.max(np.abs(ref["rho"]))),
                   **{k: tolist(x) for k, x in sol.items()})
        wall = os.path.join(case_dir, "wall_time.txt")
        rec["wall_time"] = float(open(wall).read()) if os.path.exists(wall) else None
        summary["three_d"][f"nr{nr}"] = rec
        rows.append(table_row(f"3D nr = {nr}", e)
                    + f"   (non-homogeneity {rec['non_homogeneity']:.1e}"
                    + (f", {rec['wall_time']:.0f} s)" if rec["wall_time"] else ")"))
        field_rows.append(table_row(f"3D nr = {nr}", fe, FIELD_QUANTITIES, 14))
    for ne in (2, 4, 8, 16, 32):  # 0D mesh convergence against the same reference
        z = zero_d(ne)
        e = errors(z, ref, L)
        fe = field_errors(dict(R=z["R"], fields=fields_of(z, "sigma_rr" in ref_fields)),
                          ref_fields, ref["R"], ref["P"], 1)
        summary["zero_d"][f"ne{ne}"] = dict(e, fields=fe)
        rows.append(table_row(f"0D ne = {ne}", e))
        field_rows.append(table_row(f"0D ne = {ne}", fe, FIELD_QUANTITIES, 14))
    print("max relative error over the ramp\n" + table_header() + "\n" + "\n".join(rows))
    print("strains (relative to max |E_RR|, |E_TT|), J - 1, lambda_f - 1, stresses (relative to\n"
          "their max; sigma_rr to P), max over the ramp\n"
          + table_header(FIELD_QUANTITIES, 14) + "\n" + "\n".join(field_rows))

    os.makedirs(out_dir, exist_ok=True)
    name = f"summary_ne{a.reference_ne}.json" if a.reference_ne else "summary.json"
    with open(os.path.join(out_dir, name), "w") as f:
        json.dump(summary, f, indent=1)
    if a.reference_ne:
        return
    twisting = has_twist(ref, L)
    plot(summary, ref, ref_fields, ref_name, mat, L, twisting, a.formulation, out_dir)


def plot(summary, ref, ref_fields, ref_name, mat, L, twisting, formulation, out_dir):
    """All comparisons as panels of CASE/output/comparison_<formulation>.pdf."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mmhg = C3.MMHG
    Pm = summary["P_max_mmHg"]
    levels = summary["three_d"]
    P0 = np.asarray(ref["P"]) / mmhg
    stresses = "sigma_rr" in ref_fields
    panels = (["pv", "rho"] + (["twist"] if twisting else []) + ["strain", "J", "lambda_f"]
              + (["stress", "sigma_rr"] if stresses else []) + ["conv", "conv_fields"])
    nc = 3 if len(panels) <= 9 else 4
    nrow = -(-len(panels) // nc)
    fig, axs = plt.subplots(nrow, nc, figsize=(5.2 * nc, 4.4 * nrow))
    ax = dict(zip(panels, axs.flat))
    for a in axs.flat[len(panels):]:
        a.axis("off")
    R = np.asarray(ref["R"]) * 1e3
    j = int(np.argmin(np.abs(P0 - Pm)))  # end of the ramp
    jr = lambda r: int(np.argmin(np.abs(np.array(r["P"]) / mmhg - Pm)))  # noqa: E731

    def bands(a):
        a.axhspan(5, 12, color="0.9", lw=0)
        a.axhspan(100, 120, color="0.9", lw=0)
        a.text(0.02, 8.5 / Pm, "end-diastolic", transform=a.transAxes, va="center")
        a.text(0.02, 110 / Pm, "peak systolic", transform=a.transAxes, va="center")
        a.set(ylabel="cavity pressure [mmHg]", ylim=(0, Pm))

    def profile(a, get, ls="-", mk="o", name=None, label_levels=True):
        """0D reference (line) and 3D levels (markers) through the wall at P_max."""
        a.plot(R, get(ref_fields, j), "k" + ls, lw=1.5, label=name or ref_name)
        for i, (k, r) in enumerate(levels.items()):
            a.plot(np.array(r["R"]) * 1e3, get(r["fields"], jr(r)), mk, ms=4, mfc="none",
                   c=f"C{i}", label=f"3D {k}" if label_levels else None)
        a.set(xlabel="reference radius R [mm]")
        a.legend(fontsize=8)

    def element_profile(a, key, scale=1.0, shift=0.0):
        """As profile, with the 3D element means (steps) over the nodal values."""
        a.plot(R, (np.asarray(ref_fields[key][j]) - shift) * scale, "k-", lw=1.5, label=ref_name)
        for i, (k, r) in enumerate(levels.items()):
            Rr = np.array(r["R"]) * 1e3
            y = (np.array(r["fields"][key][jr(r)]) - shift) * scale
            a.plot(Rr, y, "o", ms=3, mfc="none", c=f"C{i}", alpha=0.5)
            a.stairs(element_means(y[None, :], r["order"])[0], Rr[::r["order"]], baseline=None,
                     color=f"C{i}", lw=1.5, label=f"3D {k} (element mean; nodes: circles)")
        a.set(xlabel="reference radius R [mm]")
    row = lambda key, scale=1.0, shift=0.0: (  # noqa: E731
        lambda f, i: (np.asarray(f[key][i]) - shift) * scale)
    pmax = f"at {Pm:.0f} mmHg"

    # pressure-volume
    a = ax["pv"]
    bands(a)
    if mat["model"] == "neoHookean" and ref_name != "exact":  # compressible: incompressible limit
        Pe = np.linspace(0.0, Pm, 101)
        Ve = [math.pi * fixed_length_exact(p * mmhg, ref["R"][0], ref["R"][-1], mat["mu"]) ** 2 * L
              for p in Pe]
        a.plot(np.array(Ve) * 1e6, Pe, "k:", lw=1.5, label="incompressible, exact")
    a.plot(np.asarray(ref["V"]) * 1e6, P0, "k-", lw=2.5, label=ref_name)
    for i, (k, r) in enumerate(levels.items()):
        a.plot(np.array(r["V"]) * 1e6, np.array(r["P"]) / mmhg, "--", c=f"C{i}", label=f"3D {k}")
    a.set(xlabel="cavity volume [mL]", title="pressure-volume")
    a.legend(loc="lower right", fontsize=8)

    # radial displacement through the wall at end diastole and peak pressure
    a = ax["rho"]
    for P_show in (10.0, Pm):
        jj = int(np.argmin(np.abs(P0 - P_show)))
        a.plot(R, np.asarray(ref["rho"])[jj] * 1e3, "k-", lw=1.5)
        a.annotate(f"{P0[jj]:.0f} mmHg", (R[-1], ref["rho"][jj][-1] * 1e3), ha="right",
                   va="bottom")
        for i, (k, r) in enumerate(levels.items()):
            ji = int(np.argmin(np.abs(np.array(r["P"]) / mmhg - P_show)))
            a.plot(np.array(r["R"]) * 1e3, np.array(r["rho"][ji]) * 1e3, "o", ms=4, mfc="none",
                   c=f"C{i}", label=f"3D {k}" if P_show == Pm else None)
    a.plot([], [], "k-", label=ref_name)
    a.set(xlabel="reference radius R [mm]", ylabel="radial displacement rho [mm]",
          title="through-wall displacement")
    a.legend(fontsize=8)

    if twisting:
        a = ax["twist"]
        bands(a)
        a.plot(np.degrees(np.asarray(ref["beta"]) * L), P0, "k-", lw=2.5, label=ref_name)
        for i, (k, r) in enumerate(levels.items()):
            a.plot(np.degrees(np.array(r["beta"]) * L), np.array(r["P"]) / mmhg, "--", c=f"C{i}",
                   label=f"3D {k}")
        a.set(xlabel="twist over the length [deg]", title="twist")
        a.legend(loc="lower right", fontsize=8)

    # strains, Jacobian, fiber stretch, stresses through the wall at P_max
    a = ax["strain"]
    for n, (key, ls, mk) in enumerate((("E_RR", "-", "o"), ("E_TT", "--", "s"),
                                       ("E_ZZ", ":", "^"))):
        profile(a, row(key), ls, mk, f"{key}, {ref_name}", n == 0)
    a.set(ylabel="Green-Lagrange strain", title=f"strain (reference basis) {pmax}")
    element_profile(ax["J"], "J", shift=1.0)
    ax["J"].set(ylabel="J - 1", title=f"Jacobian {pmax}")
    ax["J"].legend(fontsize=8)
    profile(ax["lambda_f"], row("lambda_f"))
    ax["lambda_f"].set(ylabel="fiber stretch lambda_f", title=f"fiber stretch {pmax}")
    if stresses:
        a = ax["stress"]
        dev = lambda f, i: (np.asarray(f["sigma_tt"][i])  # noqa: E731
                            - np.asarray(f["sigma_rr"][i])) / 1e3
        profile(a, row("sigma_vm", 1e-3), "-", "o", f"von Mises, {ref_name}")
        profile(a, dev, "--", "s", f"sigma_tt - sigma_rr, {ref_name}", False)
        a.set(ylabel="Cauchy stress [kPa]", title=f"deviatoric stress {pmax}")
        a = ax["sigma_rr"]
        element_profile(a, "sigma_rr", 1e-3)
        a.plot(R[0], -Pm * mmhg / 1e3, "kx", ms=9, mew=2, label="-P (endo traction)")
        a.set(ylabel="sigma_rr [kPa]", title=f"radial stress {pmax}")
        a.legend(fontsize=8)

    # mesh convergence
    def convergence(a, quantities, key, slopes, title):
        if not levels:
            return
        h = np.array([r["h"] for r in levels.values()]) * 1e3
        for k, lab, mk in quantities:
            err = [r[key].get(k) for r in levels.values()]
            if None not in err and max(err) > 0.0:
                a.loglog(h, err, mk, label=lab)
        e0 = next(iter(levels.values()))[key][quantities[0][0]]
        for p, ls in zip(slopes, ("k:", "k-.")):
            a.loglog(h, e0 * (h / h[0]) ** p, ls, label=f"slope {p}")
        a.set(xlabel="element size h = (Re - Ri)/nr [mm]", ylabel="max relative error vs 0D",
              title=title)
        a.grid(True, which="both", alpha=0.3)
        a.legend(fontsize=7)
    convergence(ax["conv"], QUANTITIES, "errors", (2, 4),
                f"kinematics: 3D convergence ({formulation})")
    convergence(ax["conv_fields"], FIELD_QUANTITIES, "field_errors", (1, 2),
                f"strain, stress: 3D convergence ({formulation})")
    for n, a in enumerate(axs.flat[:len(panels)]):
        title = a.get_title()
        a.set_title("")
        a.set_title(f"({chr(97 + n)}) {title}", loc="left", fontsize=10)
    where = " (mid-length section)" if summary.get("section") == "mid-length" else ""
    fig.suptitle(f"{summary['case']}: 0D ({ref_name}) vs 3D svMultiPhysics{where}", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"comparison_{formulation}.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
