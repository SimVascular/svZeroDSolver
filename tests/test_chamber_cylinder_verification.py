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
    """Ends fixed_length: quasi-static inflation of an incompressible
    neo-Hookean tube vs. its closed-form plane-strain solution, at every load
    step."""
    Ri, Re, L = GEOM["Ri"], GEOM["Re"], GEOM["length"]
    mu = 1e4
    vv = dict(GEOM, C1=mu / 2, C2=0.0, C3=0.0, C4=0.0, C5=0.0, C6=0.0, gamma=0.0,
              num_elements=16, density=0.0, **CV.NOACT)
    T = 2000.0  # slow: no rate effect of the BCS fibre dashpot
    cfg = CV.isolated_config(vv, [0.0, T / 2, T], [0.0, 2e3, 2e3], T, 21,
                             {"pressure_all": 0.0, "volume:cyl": math.pi * Ri**2 * L},
                             formulation={"ends": "fixed_length"})
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
        vv = dict(GEOM, alpha_endo=0.0, alpha_epi=0.0, gamma=0.0, num_elements=ne,
                  density=0.0, **ho, **CV.NOACT)
        cfg = CV.isolated_config(vv, [0.0, T / 2, T], [0.0, 16e3, 16e3], T, 25,
                                 {"pressure_all": 0.0, "volume:cyl": math.pi * Ri**2 * L},
                                 formulation={"ends": "fixed_length"})
        out = CV.svzerod.simulate(cfg)
        t, V = CV.series(out, "volume:cyl")
        _, P = CV.series(out, "pressure:IN:cyl")
        ramp = t <= T / 2
        V_exact = [math.pi * CV.fixed_length_exact_ho(p, Ri, Re, vv) ** 2 * L for p in P[ramp]]
        err.append(np.abs(V[ramp] / V_exact - 1).max())
    assert err[1] < 5e-4
    assert 3.6 < err[0] / err[1] < 4.4


def test_prescribed_active_tension_matches_closed_form():
    """Prescribed active tension (svMultiPhysics UniformUnsteady): a sin^2 pulse
    as Fourier series of samples contracts the incompressible fixed-length tube
    with circumferential fibers at constant cavity pressure; the volume follows
    the closed form evaluated with an independent port of svMultiPhysics'
    Fourier interpolation, converging at second order in the element size."""
    Ri, Re, L = GEOM["Ri"], GEOM["Re"], GEOM["length"]
    ts = np.linspace(0.0, 1.0, 201)
    Ts = np.where((ts >= 0.1) & (ts <= 0.5), 2e4 * np.sin(np.pi * (ts - 0.1) / 0.4) ** 2, 0.0)
    T = CV.fourier_interpolation(ts, Ts, 50)
    ho = dict(ho_a=1e5, ho_b=0.0, ho_af=1e5, ho_bf=1.0, ho_as=0.0, ho_bs=1.0, ho_afs=0.0,
              ho_bfs=1.0, ho_k=0.0)
    P = 1e3
    err = []
    for ne in (4, 8):
        vv = dict(GEOM, alpha_endo=0.0, alpha_epi=0.0, gamma=0.0, num_elements=ne,
                  density=0.0, **ho, **CV.NOACT)
        cfg = CV.isolated_config(vv, [0.0, 1.0], [P, P], 1.0, 41,
                                 {"pressure_all": P, "volume:cyl": math.pi * Ri**2 * L},
                                 formulation={"ends": "fixed_length"},
                                 tension={"t": ts.tolist(), "T": Ts.tolist(),
                                          "fourier_modes": 50, "ramp": False})
        out = CV.svzerod.simulate(cfg)
        t, V = CV.series(out, "volume:cyl")
        m = t > 0.0
        V_exact = [math.pi * CV.fixed_length_exact_ho(P, Ri, Re, vv, float(T(x))) ** 2 * L
                   for x in t[m]]
        err.append(np.abs(V[m] / V_exact - 1).max())
    assert err[1] < 2e-4
    assert 3.5 < err[0] / err[1] < 4.5


def _input_error(formulation=None, extra=None, edit=None):
    """Message of the error that an inconsistent ChamberCylinder input raises.
    extra: flat parameters (G.cylinder_vessel); edit: function changing the
    vessel entry."""
    vv = dict(GEOM, C1=5e3, gamma=0.0, num_elements=2, density=0.0, **CV.NOACT)
    vv.update(extra or {})
    cfg = CV.isolated_config(vv, [0.0, 1.0], [0.0, 0.0], 1.0, 3,
                             {"pressure_all": 0.0, "volume:cyl": 1e-5},
                             formulation=formulation)
    if edit:
        edit(cfg["vessels"][0])
    with pytest.raises(RuntimeError) as err:
        CV.svzerod.simulate(cfg)
    return str(err.value)


def test_component_inputs_are_checked():
    """The formulation, material, active stress and activation inputs must be
    consistent, and the names of components and types are checked."""
    assert "unknown volumetric 'neo'" in _input_error({"volumetric": "neo"})
    assert "unknown component 'viscous'" in _input_error({"viscous": "none"})
    assert "selected with the material object" in _input_error({"passive": "genet23"})
    assert "bulk_modulus needs the volumetric law st91" in _input_error(
        None, {"bulk_modulus": 1e6})
    assert "st91 needs bulk_modulus > 0" in _input_error({"volumetric": "st91"})
    assert "Unknown material type 'neo'" in _input_error(
        edit=lambda v: v["material"].update(type="neo"))
    assert "Unknown parameter a4f defined in material" in _input_error(
        edit=lambda v: v["material"].update(a4f=1.0))
    assert "b4f and b4s must be positive" in _input_error(
        edit=lambda v: v.update(material={"type": "holzapfel_ogden", "a": 1e4, "a4f": 1e4}))
    # the BCS activation (nagumo) is a rate; uniform_unsteady takes none
    assert "takes no activation_function" in _input_error(
        edit=lambda v: v.update(active_stress={"type": "uniform_unsteady", "t": [0.0, 1.0],
                                               "T": [0.0, 0.0]}))
    assert "requires a rate-type" in _input_error(
        edit=lambda v: v.update(activation_function={"type": "half_cosine", "t_active": 0.2,
                                                     "t_twitch": 0.3}))


def _volume(material, active=None, activation=None):
    """Cavity volume of a short passive inflation with the given material and
    active stress law (default: the BCS law at rest). The viscosity of genet23
    keeps the system regular at the undeformed state, where the exponential
    energy has no stiffness."""
    vv = dict(GEOM, gamma=70.0, num_elements=3, density=0.0, sigma_max=0.0, k_0=0.0)
    ic = {"pressure_all": 0.0, "volume:cyl": math.pi * GEOM["Ri"] ** 2 * GEOM["length"]}
    cfg = CV.isolated_config(vv, [0.0, 0.5, 1.0], [0.0, 2e3, 2e3], 1.0, 11, ic, act=activation)
    cfg["vessels"][0]["material"] = material
    if active:
        cfg["vessels"][0]["active_stress"] = active
    return CV.series(CV.svzerod.simulate(cfg), "volume:cyl")[1]


def test_sphere_laws_in_the_cylinder():
    """The material and active stress laws are those of ChamberSphere. Its
    Mooney-Rivlin and exponential materials are the genet23 energy with
    renamed constants, and its strain-independent active stress contracts the
    cylinder."""
    mr = _volume({"type": "mooney_rivlin", "W1": 7e3, "W2": 50.0})
    assert np.array_equal(mr, _volume({"type": "genet23", "C1": 7e3, "C2": 50.0}))
    ex = _volume({"type": "exponential", "C0": 700.0, "C1": 2.0, "C2": 50.0, "C3": 4.0})
    assert np.array_equal(ex, _volume({"type": "genet23", "C3": 700.0, "C4": 2.0, "C5": 50.0,
                                       "C6": 4.0}))
    act = _volume({"type": "mooney_rivlin", "W1": 7e3, "W2": 50.0},
                  {"type": "strain_independent", "alpha_max": 30.0, "alpha_min": -30.0,
                   "sigma_max": 5e4},
                  {"type": "half_cosine", "t_active": 0.2, "t_twitch": 0.6})
    assert act.min() < 0.9 * mr[-1]


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
