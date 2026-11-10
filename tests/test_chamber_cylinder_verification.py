"""Verification of the ChamberCylinder block against independent solutions of the
same continuum model (no reference data from genet23's figures). Fast subset of
scripts/ChamberCylinder/ChamberCylinder_verification.py, which also holds the
derivations."""
import math
import os
import sys

import numpy as np
import pytest

pytest.importorskip("scipy")
sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "scripts", "ChamberCylinder"))
import ChamberCylinder_verification as CV  # noqa: E402

GEOM = {k: CV.G.GEOM[k] for k in ("Ri", "Re", "length", "alpha_endo", "alpha_epi")}
DEG = 180.0 / math.pi


def test_transmural_shear_vanishes_and_mirror_symmetry():
    """A rotation by pi about a radial axis maps body, fibres and loads onto
    themselves and (phi, eta) onto (-phi, -eta), so phi' = eta' = 0; mirroring
    the fibre angle flips the twist and nothing else."""
    r = CV.t1_symmetry(ne=3, ncycle=2)
    assert r["shear_phi_max_deg"] < 1e-9
    assert r["shear_eta_max_mm"] < 1e-9
    assert max(r["mirror_dV_rel"], r["mirror_dP_rel"], r["mirror_deps_rel"]) < 1e-8
    assert r["mirror_beta_sum_rel"] < 1e-8


@pytest.mark.parametrize("C5", [50.0, 0.0])
def test_static_inflation_extension_torsion(C5):
    """Passive static state vs. the semi-analytical solution (exact
    incompressibility, virtual work in (r_i, lambda, beta))."""
    mat = dict(CV.PASSIVE, C5=C5)
    ref = CV.static_reference(1000.0, mat, GEOM)
    blk = CV.static_block(1000.0, mat, GEOM, ne=12, Thold=10.0, dt=1e-2)
    assert ref["ok"]
    assert abs(blk["V"] / ref["V"] - 1) < 1e-3
    assert abs(blk["eps"] - (ref["q"][1] - 1)) < 5e-4
    assert abs(blk["beta"] - ref["q"][2]) * GEOM["length"] * DEG < 0.02


def test_fixed_length_inflation_matches_closed_form():
    """fixed_length = 1: quasi-static inflation of an incompressible neo-Hookean
    tube vs. its closed-form plane-strain solution, at every load step."""
    Ri, Re, L = GEOM["Ri"], GEOM["Re"], GEOM["length"]
    mu = 1e4
    vv = dict(GEOM, C1=mu / 2, C2=0.0, C3=0.0, C4=0.0, C5=0.0, C6=0.0, gamma=0.0,
              num_elements=16, density=0.0, fixed_length=1, **CV.NOACT)
    T = 2000.0  # slow: no rate effect of the BCS fibre dashpot
    cfg = CV.isolated_config(vv, [0.0, T / 2, T], [0.0, 2e3, 2e3], T, 21,
                             {"pressure_all": 0.0, "volume:cyl": math.pi * Ri**2 * L})
    out = CV.svzerod.simulate(cfg)
    t, V = CV.series(out, "volume:cyl")
    _, P = CV.series(out, "pressure:IN:cyl")
    _, eps = CV.series(out, "eps:cyl")
    ramp = t <= T / 2
    V_exact = [math.pi * CV.fixed_length_exact(p, Ri, Re, mu) ** 2 * L for p in P[ramp]]
    assert np.abs(V[ramp] / V_exact - 1).max() < 1e-4
    assert np.abs(eps).max() < 1e-15


@pytest.mark.parametrize("k", [0.0, 100.0])
def test_holzapfel_ogden_fixed_length_inflation_matches_closed_form(k):
    """Holzapfel-Ogden law (svMultiPhysics struct/LV_HolzapfelOgden_passive
    parameters in SI) with circumferential fibers and radial sheets: inflation
    of the incompressible fixed-length tube to 120 mmHg vs. its plane-strain
    closed form, converging at second order in the 0D element size."""
    Ri, Re, L = GEOM["Ri"], GEOM["Re"], GEOM["length"]
    ho = dict(ho_a=59.0, ho_b=8.023, ho_af=18472.0, ho_bf=16.026, ho_as=2481.0, ho_bs=11.12,
              ho_afs=216.0, ho_bfs=11.436, ho_k=k)
    T = 2400.0
    err = []
    for ne in (8, 16):
        vv = dict(GEOM, alpha_endo=0.0, alpha_epi=0.0, C1=0.0, C2=0.0, C3=0.0, C4=0.0, C5=0.0,
                  C6=0.0, gamma=0.0, num_elements=ne, density=0.0, fixed_length=1, **ho,
                  **CV.NOACT)
        cfg = CV.isolated_config(vv, [0.0, T / 2, T], [0.0, 16e3, 16e3], T, 25,
                                 {"pressure_all": 0.0, "volume:cyl": math.pi * Ri**2 * L})
        out = CV.svzerod.simulate(cfg)
        t, V = CV.series(out, "volume:cyl")
        _, P = CV.series(out, "pressure:IN:cyl")
        ramp = t <= T / 2
        V_exact = [math.pi * CV.fixed_length_exact_ho(p, Ri, Re, vv) ** 2 * L for p in P[ramp]]
        err.append(np.abs(V[ramp] / V_exact - 1).max())
    assert err[1] < 5e-4
    assert 3.6 < err[0] / err[1] < 4.4


def test_torsional_free_vibration():
    """Small-amplitude torsion of the free tube vs. the closed-form frequency and
    decay (twist inertia with momentum-free rigid rotation, shear stiffness,
    tissue viscosity and the fibre dashpot of the unactivated active branch)."""
    wd, s, _ = CV.twist_reference(1e4, 50.0, 70.0, 1000.0)
    t, b = CV.vibration_block("twist", dt=1e-4, T=0.15, gamma=50.0, mu=70.0)
    fit = CV.damped_fit(t, b, wd, s)
    assert abs(fit["omega"] / wd - 1) < 1e-3
    assert abs(fit["sigma"] / s - 1) < 5e-3


def test_unit_scaling_invariance():
    """The same torsion test in (mm, kPa, s) units gives the same history."""
    _, b1 = CV.vibration_block("twist", dt=1e-4, T=0.05)
    _, b2 = CV.vibration_block("twist", scale=(1e3, 1e-3), dt=1e-4, T=0.05)
    assert np.abs(b1 - b2).max() < 1e-6 * np.abs(b1).max()
