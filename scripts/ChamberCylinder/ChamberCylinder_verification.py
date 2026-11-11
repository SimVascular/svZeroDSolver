# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Verification of the ChamberCylinder block independent of the genet23 figures.

The block's kinematics (r = R + rho(R), theta = Theta + beta Z + phi(R),
z = (1 + eps) Z + eta(R)) belong to the classical inflation-extension-torsion
family of a thick-walled tube. With helical fibres in the Theta-Z plane, a
rotation by pi about a radial axis maps the body, fibres and loads onto
themselves and the solution (phi, eta) onto (-phi, -eta), so the transmural
shears phi', eta' vanish identically. With exact incompressibility the
continuum problem then reduces to three generalized coordinates
q = (r_i, lambda, beta) plus the active internal variables, which this script
solves with independent code:

  T1  symmetry: transmural shears vanish; mirrored fibres flip the twist only
  T2  static inflation-extension-torsion vs. an independent semi-analytical
      solution (virtual work in q, exact incompressibility, Gauss-Legendre)
  T3  free vibration: torsion vs. the closed-form frequency and damping, axial
      and breathing modes vs. an independent linear modal analysis, and a
      unit-rescaling check
  T4  full dynamics: an independent closed-loop re-implementation (q, BCS ODEs,
      genet23 valves and Windkessel, Radau) vs. the block's pressure, volume
      and twist over a heart beat

Only the model *definition* is shared with the block (genet23 Eqs 8-33: strain
energy, viscosity, active law, inertia, loads); kinematics, quadrature, mass,
forces and time integration are written from scratch here.

Usage (from the repository root, after building the Python module):
  PYTHONPATH=build/python python scripts/ChamberCylinder/ChamberCylinder_verification.py \
      [T1 T2 T3 T4 plot] [--out DIR]
T4 takes a few minutes. Results (JSON) and the summary figure
verification_summary.png go to DIR (default scripts/ChamberCylinder/output/).
A fast subset runs in tests/test_chamber_cylinder_verification.py.
"""
import argparse
import json
import math
import os
import sys

import numpy as np
import svzerod
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, curve_fit, root

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ChamberCylinder_genet_exact as G  # noqa: E402  (baseline parameters only)

H = 1e-30  # complex-step size


# --------------------------------------------------------------------------- #
#  Independent reduced model (exact incompressibility, phi' = eta' = 0)
# --------------------------------------------------------------------------- #
def frank_starling(ec):
    """n0(e_c): the 9-point piecewise-linear force-length curve (Caruel 2013)."""
    x = [-0.1668, -0.0073, 0.0534, 0.0969, 0.1326, 0.2016, 0.4663, 0.9187, 1.1762]
    y = [0.0, 0.5614, 0.7748, 0.8933, 0.9618, 1.0, 1.0, 0.1075, 0.0]
    return np.interp(ec, x, y)


def nagumo(t, tsys, tdias, qrs, ramp, amax, amin, T=0.8):
    """Activation rate nu(t) (an input of the model, not under test)."""
    t = t % T
    t_max = tsys + qrs
    t_plat = tsys + qrs + max(0.0, tdias - tsys - 2 * qrs)
    if t < t_max:
        nu = amax * (t - tsys) / (t_max - tsys)
    elif t < t_plat:
        nu = amax
    elif t < tdias:
        nu = amax * (1 - (t - t_plat) / (tdias - t_plat))
    else:
        nu = amin * (t - tdias) / (tdias + ramp - tdias)
    return max(nu, amin)


class Cylinder:
    """Thick-walled cylinder, q = (r_i, lambda, beta), Gauss-Legendre in R."""

    def __init__(self, Ri, Re, L, a_endo, a_epi, C, gamma=0.0, rho0=0.0, nq=40,
                 act=None, mu_fib=0.0):
        x, w = np.polynomial.legendre.leggauss(nq)
        self.Ri, self.Re, self.L, self.C, self.gamma, self.rho0 = Ri, Re, L, C, gamma, rho0
        self.R = 0.5 * (Re + Ri) + 0.5 * (Re - Ri) * x
        self.wv = 2 * np.pi * L * 0.5 * (Re - Ri) * w * self.R  # dV = 2 pi L R dR
        frac = (self.R - Ri) / (Re - Ri)
        al = np.deg2rad(a_endo * (1 - frac) + a_epi * frac)  # genet23 Eq 22
        self.ca, self.sa = np.cos(al), np.sin(al)
        self.act = act
        self.mu_fib = mu_fib  # fibre dashpot of an unactivated BCS branch (k_s -> inf)
        self.nq = nq

    # kinematics -------------------------------------------------------------
    def r(self, q):
        ri, lam, _ = q
        return np.sqrt(ri ** 2 + (self.R ** 2 - self.Ri ** 2) / lam)  # J = 1

    def Ccomp(self, q):
        """Right Cauchy-Green tensor (orthonormal R, Theta, Z), as 3x3 arrays.
        F = [[r', 0, 0], [0, r/R, r beta], [0, 0, lambda]], r' = R/(lambda r)."""
        ri, lam, beta = q
        r = self.r(q)
        dr = self.R / (lam * r)
        z = np.zeros_like(r)
        F = np.array([[dr, z, z], [z, r / self.R, r * beta], [z, z, lam + z]])
        return np.einsum("kin,kjn->ijn", F, F)

    def invariants(self, q):
        Cm = self.Ccomp(q)
        I1 = Cm[0, 0] + Cm[1, 1] + Cm[2, 2]
        I2 = 0.5 * (I1 ** 2 - np.einsum("ijn,ijn->n", Cm, Cm))
        I4 = self.ca ** 2 * Cm[1, 1] + 2 * self.ca * self.sa * Cm[1, 2] + self.sa ** 2 * Cm[2, 2]
        return Cm, I1, I2, I4

    def U(self, q):
        """Strain energy, genet23 Eq 25 (J = 1, so reduced = full invariants)."""
        c = self.C
        _, I1, I2, I4 = self.invariants(q)
        W = (c["C1"] * (I1 - 3) + c["C2"] * (I2 - 3)
             + c["C3"] * (np.exp(c["C4"] * (I1 - 3) ** 2) - 1)
             + c["C5"] * (np.exp(c["C6"] * (I4 - 1) ** 2) - 1))
        return np.sum(self.wv * W)

    def V(self, q):
        return np.pi * q[0] ** 2 * q[1] * self.L

    @staticmethod
    def cgrad(f, q):
        """Complex-step gradient of a scalar or array function of q."""
        out = []
        for k in range(3):
            qc = np.array(q, dtype=complex)
            qc[k] += 1j * H
            out.append(np.imag(f(qc)) / H)
        return np.array(out)

    # inertia (rigid rotation/translation eliminated: momentum-free, Z' = Z - L/2)
    def M(self, q):
        ri, lam, _ = q
        r = self.r(q)
        drq = np.array([ri / r, -(self.R ** 2 - self.Ri ** 2) / (2 * lam ** 2 * r), 0 * r])
        Mm = self.rho0 * np.einsum("in,jn,n->ij", drq, drq, self.wv)
        L2 = self.L ** 2 / 12.0
        Mm[2, 2] = Mm[2, 2] + self.rho0 * L2 * np.sum(self.wv * r ** 2)  # twist
        Mm[1, 1] = Mm[1, 1] + self.rho0 * L2 * np.sum(self.wv)           # axial
        return Mm

    def conv(self, q, qd):
        """Lagrange convective terms dM/dt qd - 1/2 d(qd^T M qd)/dq."""
        dM = self.cgrad(self.M, q)  # dM[k] = dM/dq_k
        return np.einsum("kij,k,j->i", dM, qd, qd) - 0.5 * np.einsum("kij,i,j->k", dM, qd, qd)

    # generalized forces -------------------------------------------------------
    def forces(self, q, qd, P, ec=None, tauc=None):
        dC = self.cgrad(lambda x: self.Ccomp(x), q)            # (3, 3, 3, nq)
        Edot = 0.5 * np.einsum("kijn,k->ijn", dC, qd)
        Qv = self.gamma * np.einsum("ijn,kijn,n->k", Edot, 0.5 * dC, self.wv)  # Eq 28
        Qe = self.cgrad(self.U, q)
        Qp = P * self.cgrad(self.V, q)
        Qa = np.zeros(3)
        if ec is not None:  # active fibre stress, Eqs 29-30
            _, _, _, I4 = self.invariants(q)
            lf = np.sqrt(I4)
            s1d = self.act["k_s"] * (lf - 1 - ec) / lf
            dI4 = self.cgrad(lambda x: self.invariants(x)[3], q)
            Qa = np.einsum("n,kn,n->k", s1d, 0.5 * dI4, self.wv)
        elif self.mu_fib:  # tau_c = k_c = 0, e_c = e_fib: sigma_1D = mu e_fib_dot / sqrt(I4)
            _, _, _, I4 = self.invariants(q)
            dI4 = self.cgrad(lambda x: self.invariants(x)[3], q)
            s1d = self.mu_fib * np.einsum("kn,k->n", dI4, qd) / (2 * I4)
            Qa = np.einsum("n,kn,n->k", s1d, 0.5 * dI4, self.wv)
        return Qp - Qe - Qv - Qa

    def bcs_rhs(self, q, ec, tauc, kc, nu):
        """BCS law, Eqs 32-33 (|e_c dot| smoothed by fv_reg as in the block)."""
        a = self.act
        _, _, _, I4 = self.invariants(q)
        efib = np.sqrt(I4) - 1
        ecd = (a["k_s"] * (efib - ec) - tauc) / a["mu"]
        decay = abs(nu) + a["alpha"] * np.sqrt(ecd ** 2 + a["fv_reg"])
        n0 = frank_starling(ec)
        taud = -decay * tauc + n0 * a["sigma0"] * max(nu, 0.0) + kc * ecd
        kcd = -decay * kc + n0 * a["k0"] * max(nu, 0.0)
        return ecd, taud, kcd


# --------------------------------------------------------------------------- #
#  svZeroD helpers
# --------------------------------------------------------------------------- #
def series(out, name):
    d = out[out["name"] == name].sort_values("time")
    return d["time"].to_numpy(), d["y"].to_numpy()


def isolated_config(vv, P_t, P_v, T, npts, ic, act=None, integrator="stiff"):
    """Cylinder alone: prescribed cavity pressure at the inlet, closed outlet."""
    act = act or {"type": "nagumo", "tsys": 0.1, "tdias": 0.4, "qrs": 0.08,
                  "ramp": 0.01, "alpha_max": 0.0, "alpha_min": 0.0}
    return {
        "boundary_conditions": [
            {"bc_name": "IN", "bc_type": "PRESSURE", "bc_values": {"P": P_v, "t": P_t}},
            {"bc_name": "OUT", "bc_type": "FLOW", "bc_values": {"Q": [0.0, 0.0], "t": [0.0, T]}}],
        "simulation_parameters": {
            "number_of_cardiac_cycles": 1, "number_of_time_pts_per_cardiac_cycle": npts,
            "cardiac_period": T, "steady_initial": False, "output_variable_based": True,
            "output_all_cycles": True,
            "absolute_tolerance": 1e-9, "maximum_nonlinear_iterations": 50,
            "integrator": integrator, "rho_infty": 0.5},
        "initial_condition": ic,
        "vessels": [{"boundary_conditions": {"inlet": "IN", "outlet": "OUT"}, "vessel_id": 0,
                     "vessel_length": 1.0, "vessel_name": "cyl",
                     "zero_d_element_type": "ChamberCylinder", "zero_d_element_values": vv,
                     "activation_function": act}],
        "junctions": [], "valves": []}


PASSIVE = dict(C1=7.0, C2=0.0, C3=700.0, C4=2.0, C5=50.0, C6=4.0)
NOACT = dict(sigma_max=0.0, k_0=0.0)


# --------------------------------------------------------------------------- #
#  T1 symmetry
# --------------------------------------------------------------------------- #
def t1_symmetry(ne=6, ncycle=4):
    res = {}
    runs = {}
    for sgn in (1, -1):
        cfg = G.build(ne=ne, ncycle=ncycle, P_vs=G.P_VS)
        vv = cfg["vessels"][1]["zero_d_element_values"]
        vv["alpha_endo"], vv["alpha_epi"] = 60.0 * sgn, -60.0 * sgn
        cfg["simulation_parameters"]["absolute_tolerance"] = 1e-8
        out = svzerod.simulate(cfg)
        _, V = series(out, "volume:ventricle")
        _, P = series(out, "pressure:ventricle:aortic")
        _, b = series(out, "beta:ventricle")
        _, e = series(out, "eps:ventricle")
        phi = np.array([series(out, f"phi_{a}:ventricle")[1] for a in range(ne + 1)])
        eta = np.array([series(out, f"eta_{a}:ventricle")[1] for a in range(ne + 1)])
        runs[sgn] = dict(V=V, P=P, b=b, e=e)
        if sgn == 1:
            L = vv["length"]
            res["shear_phi_max_deg"] = float(np.abs(phi - phi[0]).max() * 180 / np.pi)
            res["shear_eta_max_mm"] = float(np.abs(eta - eta[0]).max() * 1e3)
            res["rigid_phi_max_deg"] = float(np.abs(phi[0]).max() * 180 / np.pi)
            res["twist_max_deg"] = float(np.abs(b).max() * L * 180 / np.pi)
    a, m = runs[1], runs[-1]
    res["mirror_dV_rel"] = float(np.abs(a["V"] - m["V"]).max() / np.abs(a["V"]).max())
    res["mirror_dP_rel"] = float(np.abs(a["P"] - m["P"]).max() / np.abs(a["P"]).max())
    res["mirror_deps_rel"] = float(np.abs(a["e"] - m["e"]).max() / np.abs(a["e"]).max())
    res["mirror_beta_sum_rel"] = float(np.abs(a["b"] + m["b"]).max() / np.abs(a["b"]).max())
    res["pass"] = (res["shear_phi_max_deg"] < 1e-9 and res["shear_eta_max_mm"] < 1e-9
                   and max(res["mirror_dV_rel"], res["mirror_dP_rel"], res["mirror_deps_rel"],
                           res["mirror_beta_sum_rel"]) < 1e-6)
    return res


# --------------------------------------------------------------------------- #
#  T2 static inflation-extension-torsion
# --------------------------------------------------------------------------- #
def static_reference(P, mat, geom, nq=64, nstep=50):
    """Virtual work dU/dq = P dV/dq, solved with pressure continuation."""
    cyl = Cylinder(geom["Ri"], geom["Re"], geom["length"], geom["alpha_endo"],
                   geom["alpha_epi"], mat, nq=nq)
    sc = np.array([geom["Ri"], 1.0, 1.0 / geom["length"]])  # unknowns scaled to O(1)
    x, ok = np.array([1.0, 1.0, 0.0]), True
    for Pk in np.linspace(P / nstep, P, nstep):
        f = lambda x: cyl.forces(x * sc, np.zeros(3), Pk) * sc / (Pk * cyl.V(x * sc))
        sol = root(f, x, method="lm", tol=1e-15)
        ok = ok and float(np.abs(f(sol.x)).max()) < 1e-10
        x = sol.x
    q = x * sc
    return dict(q=q, V=cyl.V(q), ok=ok, res=float(np.abs(f(x)).max()))


def fixed_length_exact(P, Ri, Re, mu):
    """Inner radius of a pressurized fixed-length incompressible neo-Hookean tube.

    Plane strain (lambda_z = 1), W = mu/2 (I1 - 3): with lambda = r/R and
    r^2 - R^2 = r_i^2 - R_i^2, equilibrium integrates to
    P = mu [ln(l_i/l_e) - 1/(2 l_i^2) + 1/(2 l_e^2)],
    l_e^2 = 1 + (l_i^2 - 1) R_i^2/R_e^2. Mooney-Rivlin gives the same with
    mu = 2 (C1 + C2), since I1 = I2 in plane strain. Returns r_i.
    """
    def residual(li):
        le = math.sqrt(1.0 + (li * li - 1.0) * Ri * Ri / (Re * Re))
        return mu * (math.log(li / le) - 0.5 / li ** 2 + 0.5 / le ** 2) - P
    if P == 0.0:
        return Ri
    hi = 1.0 + P / mu
    while residual(hi) < 0.0:
        hi *= 2.0
    return Ri * brentq(residual, 1.0, hi, xtol=1e-15, rtol=1e-15)


def ho_family(E, a, b, k):
    """dW/dI of a Holzapfel-Ogden fiber or sheet term (as in ChamberCylinder.cpp)."""
    if a == 0.0:
        return 0.0
    o = 1.0 / (math.exp(k * E) + 1.0)
    rexp = math.exp(b * E * E)
    return a * ((1.0 - o) * E * rexp + 0.5 * k * (o - o * o) / b * (rexp - 1.0))


def fourier_interpolation(t, v, n_modes, ramp=False):
    """svMultiPhysics FourierInterpolation of samples (t_i, v_i) as a function
    of time: linear trend plus the Fourier series of the piecewise-linear
    interpolant minus the trend, periodic in t_last - t_0 (independent Python
    port of FourierInterpolation.cpp)."""
    t, v = np.asarray(t, float), np.asarray(v, float)
    T = t[-1] - t[0]
    slope = (v[-1] - v[0]) / T
    ts = t - t[0]
    vs = v - (v[0] + slope * ts)
    n = np.arange(1 if ramp else n_modes)
    re, im = np.zeros(len(n)), np.zeros(len(n))
    re[0] = np.sum(0.5 * np.diff(ts) * (vs[1:] + vs[:-1])) / T
    s = np.diff(vs) / np.diff(ts)
    for k in n[1:]:
        ko, kn = 2 * np.pi * k * ts[:-1] / T, 2 * np.pi * k * ts[1:] / T
        re[k] = 0.5 * T / (np.pi * k) ** 2 * np.sum(s * (np.cos(kn) - np.cos(ko)))
        im[k] = -0.5 * T / (np.pi * k) ** 2 * np.sum(s * (np.sin(kn) - np.sin(ko)))

    def f(time):
        x = np.clip(time - t[0], 0.0, T) if ramp else np.fmod(time - t[0], T)
        val = v[0] + x * slope
        if not ramp:
            K = 2 * np.pi * np.multiply.outer(x, n) / T
            val = val + np.sum(re * np.cos(K) - im * np.sin(K), axis=-1)
        return val
    return f


def fixed_length_exact_ho(P, Ri, Re, v, tension=0.0):
    """Inner radius of a pressurized fixed-length incompressible tube with the
    Holzapfel-Ogden law (ho_* values of v), circumferential fibers and radial
    sheets (alpha = 0), and a prescribed active fiber tension T (2nd
    Piola-Kirchhoff stress T e_F (x) e_F).

    Plane strain with lambda = r/R, r^2 - R^2 = r_i^2 - R_i^2:
    P = int_Ri^Re (sigma_tt - sigma_rr) R/r^2 dR with
    sigma_tt - sigma_rr = 2 W_1 (l^2 - l^-2) + 2 W_4f l^2 - 2 W_4s l^-2
    + T l^2, the fiber-sheet term vanishing (I8 = 0). Returns r_i.
    """
    h = {k: float(v.get("ho_" + k, 0.0)) for k in ("a", "b", "af", "bf", "as", "bs", "k")}

    def dsigma(R, ri):
        l2 = (R * R + ri * ri - Ri * Ri) / (R * R)
        I1 = l2 + 1.0 / l2 + 1.0
        w1 = 0.5 * h["a"] * math.exp(h["b"] * (I1 - 3.0))
        wf = ho_family(l2 - 1.0, h["af"], h["bf"], h["k"])
        ws = ho_family(1.0 / l2 - 1.0, h["as"], h["bs"], h["k"])
        return (2.0 * w1 * (l2 - 1.0 / l2) + 2.0 * wf * l2 - 2.0 * ws / l2
                + tension * l2) / l2 / R

    def residual(li):
        ri = li * Ri
        return quad(dsigma, Ri, Re, args=(ri,), epsabs=1e-10, epsrel=1e-12, limit=200)[0] - P
    if P == 0.0 and tension == 0.0:
        return Ri
    lo, hi = 0.99, 1.01  # residual increases with r_i
    while residual(hi) < 0.0:
        hi = 1.0 + 2.0 * (hi - 1.0)
    while residual(lo) > 0.0:
        lo = 1.0 - 2.0 * (1.0 - lo)
    return Ri * brentq(residual, lo, hi, xtol=1e-15, rtol=1e-15)


def static_block(P, mat, geom, ne, Tramp=5.0, Thold=25.0, dt=5e-3):
    vv = dict(geom)
    vv.update(mat)
    vv.update(NOACT)
    vv.update(dict(gamma=1.0, num_elements=ne, density=0.0))
    T = Tramp + Thold
    cfg = isolated_config(vv, [0.0, Tramp, T], [0.0, P, P], T, int(round(T / dt)) + 1,
                          {"pressure_all": 0.0,
                           "volume:cyl": math.pi * geom["Ri"] ** 2 * geom["length"]})
    out = svzerod.simulate(cfg)
    _, V = series(out, "volume:cyl")
    _, b = series(out, "beta:cyl")
    _, e = series(out, "eps:cyl")
    _, r0 = series(out, "rho_0:cyl")
    # the pressure input is periodic: read the state 1 s before it wraps to P(0)
    i, j = -int(round(1.0 / dt)), -int(round(5.0 / dt))
    return dict(V=V[i], eps=e[i], beta=b[i], ri=geom["Ri"] + r0[i],
                settle=float(abs(V[i] - V[j]) / V[i]))


def t2_static(Ps=(500.0, 1000.0, 2000.0), nes=(3, 6, 12, 24)):
    geom = {k: G.GEOM[k] for k in ("Ri", "Re", "length", "alpha_endo", "alpha_epi")}
    cases = {"fibres +-60": (PASSIVE, geom),
             "isotropic (C5=0)": (dict(PASSIVE, C5=0.0), geom)}
    res = {"rows": []}
    worst = 0.0
    for cname, (mat, gm) in cases.items():
        for P in Ps:
            ref = static_reference(P, mat, gm)
            L = gm["length"]
            for ne in nes:
                blk = static_block(P, mat, gm, ne)
                row = dict(case=cname, P=P, ne=ne,
                           V_ref_mL=ref["V"] * 1e6, V_mL=blk["V"] * 1e6,
                           eps_ref=ref["q"][1] - 1, eps=blk["eps"],
                           tw_ref_deg=ref["q"][2] * L * 180 / np.pi,
                           tw_deg=blk["beta"] * L * 180 / np.pi,
                           ref_ok=ref["ok"], settle=blk["settle"])
                row["dV_rel"] = abs(row["V_mL"] - row["V_ref_mL"]) / row["V_ref_mL"]
                row["dtw_deg"] = abs(row["tw_deg"] - row["tw_ref_deg"])
                res["rows"].append(row)
                if ne == max(nes):
                    worst = max(worst, row["dV_rel"], row["dtw_deg"] / 10)
    # observed order (fibres, P = 1 kPa)
    rr = [r for r in res["rows"] if r["case"] == "fibres +-60" and r["P"] == 1000.0]
    eV = [r["dV_rel"] for r in rr]
    res["order_V"] = [float(np.log2(eV[i] / eV[i + 1])) for i in range(len(eV) - 1)]
    res["pass"] = worst < 1e-3
    return res


# --------------------------------------------------------------------------- #
#  T3 free vibration
# --------------------------------------------------------------------------- #
def damped_fit(t, y, w0, s0):
    f = lambda t, A, s, w, ph, c: A * np.exp(-s * t) * np.cos(w * t + ph) + c
    A0 = y[0] - y[-1]
    p, _ = curve_fit(f, t, y, p0=[A0, s0, w0, 0.0, y[-1]], maxfev=20000)
    return dict(omega=abs(p[2]), sigma=p[1], fit_rms=float(np.sqrt(np.mean((f(t, *p) - y) ** 2)) / abs(A0)))


def vibration_block(mode, scale=None, dt=5e-5, T=0.3, ne=6, C1=1e4, gamma=5.0, rho0=1000.0,
                    mu=70.0):
    """Isotropic neo-Hookean tube, zero cavity pressure, released from a small
    twist (mode 'twist') or a small uniaxial stretch (mode 'axial')."""
    Ri, Re, L = G.GEOM["Ri"], G.GEOM["Re"], G.GEOM["length"]
    # unit system: lengths x sl, pressures x sp, time unchanged -> mass x sp sl
    sl, sp = (1.0, 1.0) if scale is None else scale
    vv = dict(Ri=Ri * sl, Re=Re * sl, length=L * sl, alpha_endo=60.0, alpha_epi=-60.0,
              C1=C1 * sp, C2=0.0, C3=0.0, C4=1.0, C5=0.0, C6=1.0, gamma=gamma * sp,
              num_elements=ne, density=rho0 * sp / sl ** 2, sigma_max=0.0, k_0=0.0,
              k_s=1e8 * sp, mu=mu * sp)
    ic = {"pressure_all": 0.0, "volume:cyl": math.pi * (Ri * sl) ** 2 * L * sl}
    if mode == "twist":
        ic["beta:cyl"] = 0.05 / sl
    else:
        e0 = 2e-3
        ic["eps:cyl"] = e0
        Rn = np.linspace(Ri, Re, ne + 1) * sl
        for a, Ra in enumerate(Rn):
            ic[f"rho_{a}:cyl"] = Ra * ((1 + e0) ** -0.5 - 1)
        ic["volume:cyl"] = math.pi * (Ri * sl * (1 + e0) ** -0.5) ** 2 * L * sl * (1 + e0)
    cfg = isolated_config(vv, [0.0, T], [0.0, 0.0], T, int(round(T / dt)) + 1, ic,
                          integrator="")
    cfg["simulation_parameters"].pop("integrator")
    out = svzerod.simulate(cfg)
    t, b = series(out, "beta:cyl")
    _, e = series(out, "eps:cyl")
    return t, (b * sl if mode == "twist" else e)


def twist_reference(C1, gamma, mu, rho0, geom=None):
    """Closed-form small-amplitude torsion of the free tube (rigid rotation
    momentum-free): M = rho0 Jp L^3/12, K = 2 C1 L Jp, D = gamma L Jp/2 plus the
    fibre dashpot mu int (sin a cos a R)^2 dV (e_FF = sin a cos a R beta)."""
    g = geom or G.GEOM
    Ri, Re, L = g["Ri"], g["Re"], g["length"]
    Jp = np.pi * (Re ** 4 - Ri ** 4) / 2
    x, w = np.polynomial.legendre.leggauss(64)
    R = 0.5 * (Re + Ri) + 0.5 * (Re - Ri) * x
    a = np.deg2rad(g["alpha_endo"] + (g["alpha_epi"] - g["alpha_endo"]) * (R - Ri) / (Re - Ri))
    Dmu = mu * 2 * np.pi * L * np.sum(0.5 * (Re - Ri) * w * (np.sin(a) * np.cos(a)) ** 2 * R ** 3)
    M, K, D = rho0 * Jp * L ** 3 / 12, 2 * C1 * L * Jp, gamma * L * Jp / 2 + Dmu
    s = D / (2 * M)
    return math.sqrt(K / M - s ** 2), s, Dmu / (gamma * L * Jp / 2)


def two_mode_fit(t, y, modes):
    f = lambda t, A1, s1, w1, p1, A2, s2, w2, p2, c: (
        A1 * np.exp(-s1 * t) * np.cos(w1 * t + p1) + A2 * np.exp(-s2 * t) * np.cos(w2 * t + p2) + c)
    (w1, s1), (w2, s2) = modes
    A0 = y[0] - y[-1]
    p, _ = curve_fit(f, t, y, p0=[A0, s1, w1, 0, 0.1 * A0, s2, w2, 0, y[-1]], maxfev=50000)
    return p, float(np.sqrt(np.mean((f(t, *p) - y) ** 2)) / abs(A0))


def t3_vibration(out=None):
    Ri, Re, L = G.GEOM["Ri"], G.GEOM["Re"], G.GEOM["length"]
    C1, rho0 = 1e4, 1000.0
    res = {}
    curves = {}
    # torsion, two damping splits (fibre dashpot dominant / tissue viscosity dominant)
    for gamma, mu in ((5.0, 70.0), (50.0, 70.0)):
        wd, s, ratio = twist_reference(C1, gamma, mu, rho0)
        for dt in (1e-4, 5e-5):
            t, b = vibration_block("twist", dt=dt, gamma=gamma, mu=mu)
            fit = damped_fit(t, b, wd, s)
            if dt == 5e-5:
                curves[f"g{gamma:g}"] = dict(t=t.tolist(), beta=b.tolist(), wd=wd, s=s)
            res[f"twist_g{gamma:g}_mu{mu:g}_dt{dt:g}"] = dict(
                omega=fit["omega"], omega_ref=wd, sigma=fit["sigma"], sigma_ref=s,
                Dmu_over_Dgamma=ratio, fit_rms=fit["fit_rms"],
                dw_rel=abs(fit["omega"] - wd) / wd, ds_rel=abs(fit["sigma"] - s) / s)
    # axial + breathing: independent linear modal analysis of (r_i, lambda)
    gamma, mu = 5.0, 70.0
    cyl = Cylinder(Ri, Re, L, 60.0, -60.0, dict(C1=C1, C2=0.0, C3=0.0, C4=1.0, C5=0.0, C6=1.0),
                   gamma=gamma, rho0=rho0, mu_fib=mu)
    q0 = np.array([Ri, 1.0, 0.0])
    K = np.zeros((3, 3)); D = np.zeros((3, 3))
    for j in range(3):
        dq = np.zeros(3); dq[j] = 1e-7 * (Ri if j == 0 else 1.0)
        K[:, j] = -(cyl.forces(q0 + dq, np.zeros(3), 0.0) - cyl.forces(q0 - dq, np.zeros(3), 0.0)) / (2 * dq[j])
        D[:, j] = -(cyl.forces(q0, dq, 0.0) - cyl.forces(q0, -dq, 0.0)) / (2 * dq[j])
    Mq = cyl.M(q0)
    ix = np.ix_([0, 1], [0, 1])
    A = np.block([[np.zeros((2, 2)), np.eye(2)],
                  [-np.linalg.solve(Mq[ix], K[ix]), -np.linalg.solve(Mq[ix], D[ix])]])
    ev = np.linalg.eigvals(A)
    ev = sorted(ev[np.imag(ev) > 0], key=lambda z: np.imag(z))
    modes = [(float(np.imag(z)), float(-np.real(z))) for z in ev]
    t, e = vibration_block("axial", dt=5e-5, gamma=gamma, mu=mu)
    p, rms = two_mode_fit(t, e, modes)
    got = sorted([(abs(p[2]), p[1]), (abs(p[6]), p[5])])
    res["axial_breathing"] = dict(
        modes_ref=modes, modes_fit=[(float(w), float(s)) for w, s in got], fit_rms=rms,
        dw_rel=max(abs(g[0] - m[0]) / m[0] for g, m in zip(got, modes)),
        ds_rel=max(abs(g[1] - m[1]) / m[1] for g, m in zip(got, modes)))
    # units: the torsion test in (mm, kPa, s)
    t1, b1 = vibration_block("twist", dt=1e-4)
    t2, b2 = vibration_block("twist", scale=(1e3, 1e-3), dt=1e-4)
    res["units_twist_maxdiff_rel"] = float(np.abs(b1 - b2).max() / np.abs(b1).max())
    tw = [v for k, v in res.items() if k.startswith("twist") and k.endswith("5e-05")]
    worst = max([max(v["dw_rel"], v["ds_rel"]) for v in tw]
                + [res["axial_breathing"]["dw_rel"], res["axial_breathing"]["ds_rel"]])
    res["pass"] = bool(worst < 1e-2 and res["units_twist_maxdiff_rel"] < 1e-6)
    if out:
        json.dump(curves, open(os.path.join(out, "T3_curves.json"), "w"))
    return res


# --------------------------------------------------------------------------- #
#  T4 independent dynamic re-implementation
# --------------------------------------------------------------------------- #
def block_beat(ne, npts=801, ncycle=12):
    cfg = G.build(ne=ne, ncycle=ncycle, P_vs=G.P_VS)
    cfg["simulation_parameters"]["number_of_time_pts_per_cardiac_cycle"] = npts
    cfg["simulation_parameters"]["absolute_tolerance"] = 1e-7
    vv = cfg["vessels"][1]["zero_d_element_values"]
    out = svzerod.simulate(cfg)
    t, V = series(out, "volume:ventricle")
    m = t >= t.max() - 0.8 - 1e-9
    d = dict(t=t[m] - t[m][0], t0=t[m][0], V=V[m])
    for n, k in (("pressure:ventricle:aortic", "P"), ("beta:ventricle", "beta"),
                 ("pressure:aortic_root:J0", "Par"), ("pressure:arterial:J1", "Pd"),
                 ("eps:ventricle", "eps"), ("rho_0:ventricle", "rho0"),
                 ("vrho_0:ventricle", "vrho0"), ("veps:ventricle", "veps"),
                 ("vbeta:ventricle", "vbeta")):
        d[k] = series(out, n)[1][m]
    nqb = 3 * ne
    Rq = []
    h = (vv["Re"] - vv["Ri"]) / ne
    for e_ in range(ne):
        for g in (-math.sqrt(0.6), 0.0, math.sqrt(0.6)):
            Rq.append(vv["Ri"] + (e_ + 0.5 + 0.5 * g) * h)
    d["Rq"] = np.array(Rq)
    for v in ("ec", "tauc", "kc"):
        d[v] = np.array([series(out, f"{v}_{q}:ventricle")[1][m] for q in range(nqb)])
    d["vv"] = vv
    d["act"] = cfg["vessels"][1]["activation_function"]
    return d


def cavity_pressure(Vdot, Pat, Par):
    """Ventricular pressure from the flow balance Vdot = Q_mitral - Q_aortic with
    switching valves (R = Rmin forward, Rmax backward; genet23 Eq 35). The
    balance is continuous and strictly decreasing in P, so exactly one valve
    state is consistent (up to ties at the switching points)."""
    best = None
    for m_open in (True, False):
        for a_open in (True, False):
            Rm = G.K_at_inv if m_open else G.K_iso_inv
            Ra = G.K_ar_inv if a_open else 1.0e12
            P = (Pat / Rm + Par / Ra - Vdot) / (1 / Rm + 1 / Ra)
            viol = max(0.0, (P - Pat) if m_open else (Pat - P)) + max(0.0, (Par - P) if a_open else (P - Par))
            if best is None or viol < best[0]:
                best = (viol, P)
    return best[1]


def independent_beat(blk, nq=24, rtol=1e-7):
    """Independent closed loop: cylinder (q, BCS) + valves + C-R-C-R Windkessel
    (genet23 Eqs 35-36), with the block's atrial pressure and activation inputs,
    started from the block's state at the beginning of the beat."""
    vv = blk["vv"]
    mat = {k: vv[k] for k in ("C1", "C2", "C3", "C4", "C5", "C6")}
    act = dict(k_s=vv["k_s"], mu=vv["mu"], alpha=vv["bcs_alpha"], fv_reg=vv["fv_reg"],
               sigma0=vv["sigma_max"], k0=vv["k_0"])
    cyl = Cylinder(vv["Ri"], vv["Re"], vv["length"], vv["alpha_endo"], vv["alpha_epi"], mat,
                   gamma=vv["gamma"], rho0=vv["density"], nq=nq, act=act)
    a = blk["act"]
    nu = lambda t: nagumo(t, a["tsys"], a["tdias"], a["qrs"], a["ramp"], a["alpha_max"], a["alpha_min"])
    tat, Pat_v = G.atrial_pressure(0.8)
    Pat = lambda t: float(np.interp((blk["t0"] + t) % 0.8, tat, Pat_v))
    n = cyl.nq

    def unpack(y):
        return y[:3], y[3:6], y[6], y[7], y[8:8 + n], y[8 + n:8 + 2 * n], y[8 + 2 * n:]

    def state(t, y):
        q, qd, Par, Pd, ec, tauc, kc = unpack(y)
        Vdot = float(cyl.cgrad(cyl.V, q) @ qd)
        return cavity_pressure(Vdot, Pat(t), Par)

    def rhs(t, y):
        q, qd, Par, Pd, ec, tauc, kc = unpack(y)
        Pv = state(t, y)
        Q = cyl.forces(q, qd, Pv, ec, tauc) - cyl.conv(q, qd)
        qdd = np.linalg.solve(cyl.M(q), Q)
        Ra = G.K_ar_inv if Par < Pv else 1.0e12
        Qa = (Pv - Par) / Ra
        dPar = (Qa - (Par - Pd) / G.Rp) / G.Car
        dPd = ((Par - Pd) / G.Rp - (Pd - G.P_VS) / G.Rd) / G.Cd
        ecd, taud, kcd = cyl.bcs_rhs(q, ec, tauc, kc, nu(blk["t0"] + t))
        return np.concatenate([qd, qdd, [dPar, dPd], ecd, taud, kcd])

    q0 = np.array([vv["Ri"] + blk["rho0"][0], 1 + blk["eps"][0], blk["beta"][0]])
    qd0 = np.array([blk["vrho0"][0], blk["veps"][0], blk["vbeta"][0]])
    I4 = cyl.invariants(q0)[3]
    tau0 = np.interp(cyl.R, blk["Rq"], blk["tauc"][:, 0])
    kc0 = np.interp(cyl.R, blk["Rq"], blk["kc"][:, 0])
    ec0 = np.sqrt(I4) - 1 - tau0 / act["k_s"]
    y0 = np.concatenate([q0, qd0, [blk["Par"][0], blk["Pd"][0]], ec0, tau0, kc0])
    m = 8 + 3 * n
    sp = np.zeros((m, m))
    sp[:8, :] = 1
    sp[8:, :8] = 1
    for j in range(n):
        for b in range(3):
            for c in range(3):
                sp[8 + b * n + j, 8 + c * n + j] = 1
    atol = np.concatenate([[1e-9, 1e-8, 1e-6, 1e-7, 1e-6, 1e-4, 1e-2, 1e-2],
                           np.full(n, 1e-8), np.full(n, 1e-2), np.full(n, 1e-2)])
    sol = solve_ivp(rhs, (0.0, 0.8), y0, method="Radau", t_eval=blk["t"], rtol=rtol,
                    atol=atol, jac_sparsity=sp, max_step=1e-3)
    q = sol.y[:3]
    V = np.pi * q[0] ** 2 * q[1] * vv["length"]
    P = np.array([state(t, sol.y[:, i]) for i, t in enumerate(sol.t)])
    return dict(t=sol.t, V=V, P=P, eps=q[1] - 1, beta=q[2], ok=bool(sol.success),
                nfev=int(sol.nfev))


def _t4_job(spec):
    ne, dt = spec
    blk = block_beat(ne, npts=int(round(0.8 / dt)) + 1)
    return ne, dt, blk, independent_beat(blk)


def t4_dynamics(specs=((6, 1e-3), (6, 2.5e-4), (12, 1e-3)), out=None):
    from multiprocessing import Pool
    res = {"rows": []}
    L = G.GEOM["length"]
    curves = {}
    with Pool(len(specs)) as pool:
        jobs = pool.map(_t4_job, specs)
    for ne, dt, blk, ind in jobs:
        if True:
            tw_b = blk["beta"] * L * 180 / np.pi
            tw_i = ind["beta"] * L * 180 / np.pi
            row = dict(ne=ne, dt_ms=dt * 1e3, ok=ind["ok"],
                       dV_max_mL=float(np.abs(blk["V"] - ind["V"]).max() * 1e6),
                       dP_max_kPa=float(np.abs(blk["P"] - ind["P"]).max() * 1e-3),
                       deps_max=float(np.abs(blk["eps"] - ind["eps"]).max()),
                       dtw_max_deg=float(np.abs(tw_b - tw_i).max()),
                       tw550_block=float(np.interp(0.55, blk["t"], tw_b)),
                       tw550_indep=float(np.interp(0.55, ind["t"], tw_i)),
                       V_range_mL=float(np.ptp(blk["V"]) * 1e6), tw_range_deg=float(np.ptp(tw_b)))
            res["rows"].append(row)
            curves[f"ne{ne}_dt{dt*1e3:g}"] = dict(t=blk["t"].tolist(), V_b=(blk["V"] * 1e6).tolist(),
                                                  V_i=(ind["V"] * 1e6).tolist(), tw_b=tw_b.tolist(),
                                                  tw_i=tw_i.tolist(), P=(blk["P"] * 1e-3).tolist(),
                                                  P_i=(ind["P"] * 1e-3).tolist())
    best = min(res["rows"], key=lambda r: r["dtw_max_deg"])
    res["pass"] = all(r["ok"] for r in res["rows"]) and best["dV_max_mL"] < 0.5 and best["dtw_max_deg"] < 0.2
    if out:
        json.dump(curves, open(os.path.join(out, "T4_curves.json"), "w"))
    return res


# --------------------------------------------------------------------------- #
#  summary figure
# --------------------------------------------------------------------------- #
def plot_summary(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.4))
    L = G.GEOM["length"]
    # T2: static P-V and twist, reference curve + block points
    T2 = json.load(open(os.path.join(out, "T2.json")))
    geom = {k: G.GEOM[k] for k in ("Ri", "Re", "length", "alpha_endo", "alpha_epi")}
    Ps = np.linspace(100.0, 2000.0, 20)
    for cname, mat, col in (("fibres +-60", PASSIVE, "C0"), ("isotropic (C5=0)", dict(PASSIVE, C5=0.0), "C2")):
        refs = [static_reference(P, mat, geom) for P in Ps]
        ax[0, 0].plot([r["V"] * 1e6 for r in refs], Ps / 1e3, color=col, label=f"semi-analytical, {cname}")
        if cname.startswith("fib"):
            ax[0, 1].plot(Ps / 1e3, [r["q"][2] * L * 180 / np.pi for r in refs], color=col, label="semi-analytical")
        for ne, mk in ((3, "s"), (24, "o")):
            rows = [r for r in T2["rows"] if r["case"] == cname and r["ne"] == ne]
            ax[0, 0].plot([r["V_mL"] for r in rows], [r["P"] / 1e3 for r in rows], mk, mfc="none", color=col,
                          label=f"block, {ne} elements" if cname.startswith("fib") else None)
            if cname.startswith("fib"):
                ax[0, 1].plot([r["P"] / 1e3 for r in rows], [r["tw_deg"] for r in rows], mk, mfc="none",
                              color=col, label=f"block, {ne} elements")
    ax[0, 0].set_xlabel("cavity volume [mL]"); ax[0, 0].set_ylabel("cavity pressure [kPa]")
    ax[0, 0].set_title("T2 static inflation (passive)"); ax[0, 0].legend(fontsize=7)
    ax[0, 1].set_xlabel("cavity pressure [kPa]"); ax[0, 1].set_ylabel("torsion [deg]")
    ax[0, 1].set_title("T2 static twist (fibres +-60)"); ax[0, 1].legend(fontsize=7)
    a = ax[0, 2]
    for cname, col in (("fibres +-60", "C0"), ("isotropic (C5=0)", "C2")):
        rows = [r for r in T2["rows"] if r["case"] == cname and r["P"] == 1000.0]
        a.loglog([r["ne"] for r in rows], [r["dV_rel"] for r in rows], "o-", color=col, label=f"volume, {cname}")
        if cname.startswith("fib"):
            a.loglog([r["ne"] for r in rows], [r["dtw_deg"] / 10 for r in rows], "s--", color=col,
                     label="twist error / 10 deg")
    ne = np.array([3, 24]); a.loglog(ne, 6e-3 * (ne / 3.0) ** -2, "k:", label="slope 2")
    from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
    a.xaxis.set_major_locator(FixedLocator([3, 6, 12, 24]))
    a.xaxis.set_major_formatter(FixedFormatter(["3", "6", "12", "24"]))
    a.xaxis.set_minor_locator(NullLocator())
    a.set_xlabel("elements through the wall"); a.set_ylabel("error vs semi-analytical (1 kPa)")
    a.set_title("T2 convergence"); a.legend(fontsize=7)
    # T3: torsional free vibration
    a = ax[1, 0]
    C3 = json.load(open(os.path.join(out, "T3_curves.json")))
    for k, col in (("g5", "C1"), ("g50", "C3")):
        c = C3[k]; t = np.array(c["t"]); b = np.array(c["beta"])
        ref = b[0] * np.exp(-c["s"] * t) * (np.cos(c["wd"] * t) + c["s"] / c["wd"] * np.sin(c["wd"] * t))
        a.plot(t * 1e3, b * L * 180 / np.pi, color=col, lw=2.5, alpha=0.5, label=f"block, gamma={k[1:]}")
        a.plot(t * 1e3, ref * L * 180 / np.pi, "k--", lw=0.9)
    a.plot([], [], "k--", lw=0.9, label="closed form")
    a.set_xlim(0, 150); a.set_xlabel("time [ms]"); a.set_ylabel("torsion [deg]")
    a.set_title("T3 torsional free vibration"); a.legend(fontsize=7)
    # T4: independent dynamics
    C4 = json.load(open(os.path.join(out, "T4_curves.json")))
    key = sorted(C4, key=lambda k: (-int(k.split("_")[0][2:]), k))[0]
    c = C4[key]; t = np.array(c["t"]) * 1e3
    ax[1, 1].plot(t, c["V_b"], color="C0", lw=2.5, alpha=0.5, label=f"block ({key})")
    ax[1, 1].plot(t, c["V_i"], "k--", lw=1, label="independent re-implementation")
    ax[1, 1].set_xlabel("time [ms]"); ax[1, 1].set_ylabel("cavity volume [mL]")
    ax[1, 1].set_title("T4 full beat: volume"); ax[1, 1].legend(fontsize=7)
    ax[1, 2].plot(t, c["tw_b"], color="C3", lw=2.5, alpha=0.5, label=f"block ({key})")
    ax[1, 2].plot(t, c["tw_i"], "k--", lw=1, label="independent re-implementation")
    ax[1, 2].set_xlabel("time [ms]"); ax[1, 2].set_ylabel("torsion [deg]")
    ax[1, 2].set_title("T4 full beat: torsion"); ax[1, 2].legend(fontsize=7)
    for a in ax.flat:
        a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "verification_summary.png"), dpi=130)
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tests", nargs="*", default=["T1", "T2", "T3", "T4"])
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "output"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fns = dict(T1=t1_symmetry, T2=t2_static, T3=lambda: t3_vibration(out=a.out),
               T4=lambda: t4_dynamics(out=a.out), plot=lambda: plot_summary(a.out))
    allres = {}
    for k in a.tests:
        print(f"--- {k}", flush=True)
        r = fns[k]()
        allres[k] = r
        print(json.dumps(r, indent=1, default=float), flush=True)
        json.dump(r, open(os.path.join(a.out, f"{k}.json"), "w"), indent=1, default=float)
    print("SUMMARY: " + "  ".join(f"{k} {'PASS' if v['pass'] else 'FAIL'}" for k, v in allres.items()
                                  if "pass" in v))


if __name__ == "__main__":
    main()
