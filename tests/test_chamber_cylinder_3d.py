"""Generator of the 3D svMultiPhysics counterpart of a ChamberCylinder model
(scripts/ChamberCylinder/verification_3d/cylinder_3d.py): mesh geometry, the
0D-to-3D mapping and the lumen-volume integral. svMultiPhysics is not run."""
import json
import math
import os
import re
import sys

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "scripts", "ChamberCylinder", "verification_3d"))
import cylinder_3d as C3  # noqa: E402

CASES = os.path.join(HERE, "..", "scripts", "ChamberCylinder", "verification_3d")
MODEL = os.path.join(CASES, "passive_isotropic", "chamber_cylinder_0d.json")
MODEL_HO = os.path.join(CASES, "passive_anisotropic", "chamber_cylinder_0d.json")
Ri, Re, L = 0.019, 0.033, 0.057


def tet_volumes(x, c):
    a, b, d, e = (x[c[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, d - a), e - a) / 6.0


@pytest.mark.parametrize("sector", [90.0, 360.0])
@pytest.mark.parametrize("order", [1, 2])
def test_mesh_geometry(sector, order):
    m = C3.mesh_cylinder(Ri, Re, L, 2, sector, order)
    x, c = m["points"], m["cells"]
    v = tet_volumes(x, c)
    assert (v > 0).all()
    wall = math.radians(sector) / 2 * (Re**2 - Ri**2) * L
    assert abs(v.sum() / wall - 1) < 2e-2  # straight-sided corners: O(h^2)
    assert len(np.unique(c)) == len(x)
    extra = {"x0", "y0"} if sector == 90.0 else set()
    assert set(m["faces"]) == {"endo", "epi", "base", "top"} | extra
    if sector == 360.0:
        assert m["ntheta"] % 4 == 0
    for name, sign in (("endo", -1.0), ("epi", 1.0)):  # normals point out of the solid
        f = m["faces"][name][0]
        n = np.cross(x[f[:, 1]] - x[f[:, 0]], x[f[:, 2]] - x[f[:, 0]])
        g = x[f[:, :3]].mean(axis=1)
        assert (sign * (n[:, 0] * g[:, 0] + n[:, 1] * g[:, 1]) > 0).all()


def test_lumen_volume_reference_and_convergence():
    V = math.pi * Ri**2 * L
    err = {o: [abs(C3.lumen_volume(m["points"], m["faces"]["endo"][0], 90.0) / V - 1)
               for m in (C3.mesh_cylinder(Ri, Re, L, n, 90.0, o) for n in (2, 4))]
           for o in (1, 2)}
    assert err[1][0] / err[1][1] > 3.5  # faceted circle: second order
    assert err[2][1] < 1e-6             # quadratic faces on the exact circle


def test_lumen_volume_of_uniform_expansion():
    m = C3.mesh_cylinder(Ri, Re, L, 2, 90.0, 2)
    x = m["points"].copy()
    x[:, :2] *= 1.1
    V = C3.lumen_volume(x, m["faces"]["endo"][0], 90.0)
    assert abs(V / (math.pi * (1.1 * Ri) ** 2 * L) - 1) < 2e-5


def test_generate_from_0d(tmp_path):
    case = C3.generate(MODEL, str(tmp_path), nr=1)
    with open(MODEL) as f:
        vals = json.load(f)["vessels"][0]["zero_d_element_values"]
    assert case["material"]["model"] == "neoHookean"
    assert case["material"]["mu"] == 2.0 * vals["C1"]
    assert case["T_ramp"] == 1200.0
    assert abs(case["P_max"] / C3.MMHG - 120.0) < 1e-9  # 0 to 120 mmHg
    xml = (tmp_path / "solver.xml").read_text()
    mu, kappa = 2.0 * vals["C1"], vals["bulk_modulus"]
    assert case["material"]["kappa"] == kappa
    assert '<Add_equation type="struct" >' in xml  # displacement-based (default)
    nu = (3 * kappa - 2 * mu) / (2 * (3 * kappa + mu))
    assert f"<Elasticity_modulus> {2.0 * mu * (1 + nu):.17g} </Elasticity_modulus>" in xml
    assert f"<Penalty_parameter> {kappa:.17g} </Penalty_parameter>" in xml
    for face in case["faces"]:
        assert (tmp_path / "mesh" / "mesh-surfaces" / f"{face}.vtu").exists()
    load = [float(x) for x in (tmp_path / "endo_pressure.dat").read_text().split()]
    assert load == [2.0, 1.0, 0.0, 0.0, case["T_ramp"], case["P_max"]]


@pytest.mark.parametrize("case_name", ["passive_anisotropic", "passive_anisotropic_full"])
def test_generate_holzapfel_ogden(tmp_path, case_name):
    """Anisotropic case: HolzapfelOgden with fiber and sheet files, eighth model
    (90-degree sector, half length), and the same stress and directions as the
    0D block."""
    import check_materials
    model = os.path.join(CASES, case_name, "chamber_cylinder_0d.json")
    case = C3.generate(model, str(tmp_path), nr=1, half_length=True)
    with open(model) as f:
        vals = json.load(f)["vessels"][0]["zero_d_element_values"]
    assert case["material"]["model"] == "HolzapfelOgden"
    assert case["sector"] == 90.0 and case["mesh_length"] == 0.5 * vals["length"]
    xml = (tmp_path / "solver.xml").read_text()
    assert '<Constitutive_model type="HolzapfelOgden">' in xml
    for key, tag in zip(C3.HO_KEYS, C3.HO_SVMP):
        assert f"<{tag}> {vals['ho_' + key]:.17g} </{tag}>" in xml
    assert xml.count("<Fiber_direction_from_file_points>") == 2  # nodal directions
    err, fib = check_materials.check(model, str(tmp_path), n=300)
    assert err < 1e-12 and fib < 1e-12


def test_helical_fibers_need_the_full_cylinder(tmp_path):
    with open(MODEL_HO) as f:
        cfg = json.load(f)
    vals = cfg["vessels"][0]["zero_d_element_values"]
    vals["alpha_endo"], vals["alpha_epi"] = 60.0, -60.0
    model = tmp_path / "helical.json"
    model.write_text(json.dumps(cfg))
    with pytest.raises(ValueError):  # the mid-length plane would be a slip end
        C3.generate(str(model), str(tmp_path / "a"), nr=1, half_length=True)
    case = C3.generate(str(model), str(tmp_path / "a"), nr=1)
    assert case["sector"] == 360.0 and case["mesh_length"] == case["length"]
    with pytest.raises(ValueError):
        C3.generate(str(model), str(tmp_path / "b"), nr=1, sector=90.0)
    # rigid-body modes: in-plane dashpots (damping 1e-3 mu / Re * dt) on both ends, u_z = 0
    xml = (tmp_path / "a" / "solver.xml").read_text()
    damping = 1e-3 * case["material"]["mu"] / case["Re"] * case["dt"]
    assert case["end_damping"] == pytest.approx(damping)
    for name in ("base", "top"):
        bcs = re.findall(rf'<Add_BC name="{name}" >(.*?)</Add_BC>', xml, flags=re.S)
        assert sorted(re.search(r"<Type> (\w+) </Type>", b).group(1) for b in bcs) == ["Dir", "Robin"]
        robin = next(b for b in bcs if "Robin" in b)
        assert "<Stiffness> 0.0 </Stiffness>" in robin and f"<Damping> {damping:.6g} </Damping>" in robin
    assert '<Add_BC name="epi"' not in xml  # no spring on the epicardium


def test_kinematics_of_twisted_inflation():
    """The 0D kinematics are recovered from the exact Z-homogeneous field
    r = R + rho(R), theta = Theta + beta Z + c, z = Z."""
    m = C3.mesh_cylinder(Ri, Re, L, 2, 360.0, 2)
    X = m["points"]
    R, Th = np.hypot(X[:, 0], X[:, 1]), np.arctan2(X[:, 1], X[:, 0])
    rho = lambda R: 1e-3 * (1.0 + 30.0 * (Re - R))  # noqa: E731
    beta, c = -1.6, 0.01
    r, th = R + rho(R), Th + beta * X[:, 2] + c
    x = np.stack([r * np.cos(th), r * np.sin(th), X[:, 2]], axis=1)
    levels = np.rint((R - Ri) / ((Re - Ri) / 4)).astype(int)
    k = C3.kinematics(X, x, levels)
    assert abs(k["beta"] / beta - 1) < 1e-12
    assert np.allclose(k["rho"], rho(Ri + (Re - Ri) * np.arange(5) / 4), rtol=0, atol=1e-15)
    assert np.abs(k["phi"]).max() < 1e-12 and np.abs(k["eta"]).max() < 1e-15
    assert k["rho_spread"] < 1e-15


def test_kinematics_of_the_mid_section():
    """Mid-section kinematics: only the middle element layer enters, beta is the
    local twist rate at L/2 (here of a non-uniform twist), and an end layer of
    the radial displacement elsewhere is ignored."""
    m = C3.mesh_cylinder(Ri, Re, L, 2, 360.0, 2)
    X, Z = m["points"], m["points"][:, 2]
    R, Th = np.hypot(X[:, 0], X[:, 1]), np.arctan2(X[:, 1], X[:, 0])
    z0, dz = 0.5 * L, 0.5 * L / m["nz"]
    rho = 1e-3 * (1.0 + 30.0 * (Re - R)) * (1.0 + 0.5 * (np.abs(Z - z0) > dz * 1.01))
    beta, gamma = -1.6, 40.0  # theta = beta (Z - z0) + gamma (Z - z0)^2
    r, th = R + rho, Th + beta * (Z - z0) + gamma * (Z - z0) ** 2
    x = np.stack([r * np.cos(th), r * np.sin(th), Z], axis=1)
    levels = np.rint((R - Ri) / ((Re - Ri) / 4)).astype(int)
    k = C3.kinematics(X, x, levels, z0=z0, dz=dz)
    assert abs(k["beta"] / beta - 1) < 1e-10
    Rk = Ri + (Re - Ri) * np.arange(5) / 4
    assert np.allclose(k["rho"], 1e-3 * (1.0 + 30.0 * (Re - Rk)), rtol=0, atol=1e-15)
    assert k["rho_spread"] < 1e-15


@pytest.mark.parametrize("order", [1, 2])
def test_struct_fields_of_homogeneous_deformation(order):
    """A linear displacement u = A X is reproduced exactly by the (curved,
    isoparametric) elements: F = I + A, J = det F and E = (F^T F - I)/2 at
    every node; the stress equals the svMultiPhysics model at that F."""
    import check_materials
    m = C3.mesh_cylinder(Ri, Re, L, 2, 90.0, order)
    X = m["points"]
    A = np.array([[0.05, 0.01, 0.0], [-0.02, 0.03, 0.01], [0.0, 0.02, -0.01]])
    fib = C3.cylinder_fibers_at_nodes(X, Ri, Re, 30.0, -30.0)
    mat = dict(model="neoHookean", E=3e5, nu=0.45, Kpen=1e6, vol="ST91")
    sf = C3.struct_fields(X, m["cells"], X @ A.T, order, mat, fib["fiber"], fib["sheet"])
    F = np.eye(3) + A
    assert np.abs(sf["Def_grad"] - F.reshape(9)).max() < 1e-12
    assert np.abs(sf["Jacobian"] - np.linalg.det(F)).max() < 1e-12
    E = 0.5 * (F.T @ F - np.eye(3))
    assert np.abs(sf["Strain"][:, 3] - E[0, 1]).max() < 1e-12
    S = check_materials.stress_3d(F, fib["fiber"][0], fib["sheet"][0], mat)
    assert np.abs(sf["Stress"][:, 0] - S[0, 0]).max() < 1e-9 * np.abs(S).max()


def test_unsupported_0d_inputs():
    with pytest.raises(ValueError):
        C3.material({"C1": 1.0, "C3": 1.0})
    with pytest.raises(ValueError):  # HolzapfelOgden replaces the C1-C6 law
        C3.material({"C1": 1.0, "ho_a": 1.0, "ho_b": 1.0, "alpha_endo": 0, "alpha_epi": 0})
    with pytest.raises(ValueError):  # HolzapfelOgden needs its isotropic term
        C3.material({"ho_af": 1.0, "ho_bf": 1.0, "alpha_endo": 0, "alpha_epi": 0})
    with pytest.raises(ValueError):
        C3.pressure_ramp([0.0, 1.0, 2.0], [0.0, 1.0, 0.5])


MODEL_ACTIVE = os.path.join(CASES, "active_prescribed", "chamber_cylinder_0d.json")


def test_coupled_circulation():
    """The ChamberCylinder becomes a cavity junction whose extra inlet is the 3D
    flow (Q = -dV/dt); impedances of the eighth model are scaled by 8."""
    import coupled_3d as CP
    with open(MODEL_ACTIVE) as f:
        cfg = json.load(f)
    c = CP.coupling_config(cfg, 8.0)
    names = {v["vessel_name"]: v for v in c["vessels"]}
    assert "cylinder" not in names and {"cavity_in", "cavity_out", "cavity_3d"} <= set(names)
    orig = {v["vessel_name"]: v for v in cfg["vessels"]}
    for n in ("upstream_vessel", "downstream_vessel"):
        vals, ov = names[n]["zero_d_element_values"], orig[n]["zero_d_element_values"]
        assert vals["R_poiseuille"] == 8.0 * ov["R_poiseuille"]
        assert vals["C"] == ov["C"] / 8.0
    valves = {v["name"]: v["params"] for v in c["valves"]}
    assert valves["inlet_valve"]["downstream_block"] == "cavity_in"
    assert valves["outlet_valve"]["upstream_block"] == "cavity_out"
    assert c["junctions"][-1]["inlet_blocks"] == ["cavity_in", "cavity_3d"]
    blk = c["external_solver_coupling_blocks"][0]
    assert (blk["type"], blk["location"], blk["connected_block"]) == ("FLOW", "inlet", "cavity_3d")
    assert c["simulation_parameters"]["coupled_simulation"] is True
    assert "volume:cylinder" not in c["initial_condition"]


def test_generate_coupled(tmp_path):
    import coupled_3d as CP
    case = CP.generate(MODEL_ACTIVE, str(tmp_path), nr=1, cycles=1)
    xml = (tmp_path / "solver.xml").read_text()
    assert "<Model> UniformUnsteady </Model>" in xml
    assert "<svZeroDSolver_block> LV </svZeroDSolver_block>" in xml
    assert "<Coupling_type> implicit </Coupling_type>" in xml
    assert '<LS type="GMRES" >' in xml  # BiCGStab drops the coupled-BC tangent
    assert case["scale"] == 8.0 and case["n_steps"] == 400
    lines = (tmp_path / "active_tension.dat").read_text().split("\n")
    with open(MODEL_ACTIVE) as f:
        at = CP.find_cylinder(json.load(f))["active_tension"]
    assert lines[0].split() == [str(len(at["t"])), str(at["fourier_modes"])]


def test_zero_d_fields_of_homogeneous_state():
    """rho = c R, phi = eta = beta = 0: F = diag(1 + c, 1 + c, 1), so
    E_RR = E_TT = ((1 + c)^2 - 1)/2, E_ZZ = 0, J = (1 + c)^2, the fiber stretch
    of the helix angle a is sqrt((1 + c)^2 cos^2 a + sin^2 a) and the stress
    is that of the 0D material law at F."""
    import check_materials
    import run_study as RS
    with open(MODEL_HO) as f:
        v = json.load(f)["vessels"][0]["zero_d_element_values"]
    R = np.linspace(Ri, Re, 5)
    c, a = 0.05, (30.0, -30.0)
    out = RS.zero_d_fields(R, [c * R], [0 * R], [0 * R], [0.0], v, a)
    lam = 1.0 + c
    assert np.allclose(out["E_RR"], 0.5 * (lam ** 2 - 1.0), rtol=1e-12)
    assert np.allclose(out["E_TT"], 0.5 * (lam ** 2 - 1.0), rtol=1e-12)
    assert np.abs(out["E_ZZ"]).max() < 1e-15
    assert np.allclose(out["J"], lam ** 2, rtol=1e-12)
    alpha = np.radians(np.linspace(*a, 5))
    assert np.allclose(out["lambda_f"], np.sqrt(lam ** 2 * np.cos(alpha) ** 2 + np.sin(alpha) ** 2))
    F = np.diag([lam, lam, 1.0])
    vals = {k: float(v.get(k, 0.0)) for k in ("C1", "C2", "C3", "C4", "C5", "C6")}
    vals.update({k: x for k, x in v.items() if k.startswith("ho_") or k == "bulk_modulus"})
    for j in range(5):
        f = np.array([0.0, np.cos(alpha[j]), np.sin(alpha[j])])
        sig = F @ check_materials.stress_0d(F, f, np.array([1.0, 0.0, 0.0]), vals) @ F / lam ** 2
        for k, i in (("sigma_rr", 0), ("sigma_tt", 1), ("sigma_zz", 2)):
            assert abs(out[k][0, j] / sig[i, i] - 1.0) < 1e-12


def test_zero_d_fields_match_the_3d_fields():
    """The 0D strain, Jacobian and fiber stretch profiles of a Z-homogeneous
    field r = R + rho(R), theta = Theta + beta Z + phi(R), z = Z + eta(R) equal
    the level means of the 3D nodal fields of the same displacement on a P2
    mesh, to the discretization error."""
    import run_study as RS
    nr, a = 4, (30.0, -30.0)
    m = C3.mesh_cylinder(Ri, Re, L, nr, 90.0, 2)
    X = m["points"]
    R, Th = np.hypot(X[:, 0], X[:, 1]), np.arctan2(X[:, 1], X[:, 0])
    rho = lambda R: 2e-3 * (1.0 + 0.5 * ((Re - R) / (Re - Ri)) ** 2)  # noqa: E731
    phi = lambda R: 0.05 * ((R - Ri) / (Re - Ri)) ** 2  # noqa: E731
    eta = lambda R: 1e-3 * np.sin((R - Ri) / (Re - Ri))  # noqa: E731
    beta = -1.6
    r, th = R + rho(R), Th + beta * X[:, 2] + phi(R)
    x = np.stack([r * np.cos(th), r * np.sin(th), X[:, 2] + eta(R)], axis=1)
    levels = np.rint((R - Ri) / ((Re - Ri) / (2 * nr))).astype(int)
    fib = C3.cylinder_fibers_at_nodes(X, Ri, Re, *a)["fiber"]
    f3 = C3.level_means(levels, C3.kinematic_fields(X, m["cells"], x - X, 2, fib),
                        C3.KINEMATIC_PROFILES)
    Rl = Ri + (Re - Ri) * np.arange(2 * nr + 1) / (2 * nr)
    f0 = RS.zero_d_fields(Rl, [rho(Rl)], [phi(Rl)], [eta(Rl)], [beta], {}, a, stresses=False)
    for k in C3.KINEMATIC_PROFILES:
        scale = np.abs(f0[k] - (1.0 if k in ("J", "lambda_f") else 0.0)).max()
        assert np.abs(f3[k] - f0[k][0]).max() < 0.03 * scale, k


def test_element_means():
    """Simpson (P2) and trapezoid (P1) means are exact for quadratic and
    linear profiles."""
    import run_study as RS
    R = np.linspace(1.0, 2.0, 5)
    assert np.allclose(RS.element_means([R ** 2], 2), [[(1.5 ** 3 - 1) / 1.5, (8 - 1.5 ** 3) / 1.5]])
    assert np.allclose(RS.element_means([R], 1), [0.5 * (R[1:] + R[:-1])])


def test_active_cycle_metrics():
    """Elliptic PV loop V = V0 + a cos(2 pi t), P = P0 + b sin(2 pi t): stroke
    work pi a b, stroke volume 2 a, dP/dt extrema +-2 pi b, and the valve events
    of a sinusoidal flow."""
    import run_active_study as RA
    t = np.linspace(0.0, 1.0, 4001)
    V0, a, P0, b = 60.0, 10.0, 50.0, 40.0
    m = RA.cycle_metrics(t, V0 + a * np.cos(2 * np.pi * t), P0 + b * np.sin(2 * np.pi * t))
    assert abs(m["W"] / (np.pi * a * b) - 1.0) < 1e-6
    assert abs(m["SV"] - 2 * a) < 1e-9 and abs(m["EF"] - 2 * a / (V0 + a)) < 1e-9
    assert abs(m["dPdt+"] / (2 * np.pi * b) - 1.0) < 1e-5
    assert abs(m["dPdt-"] / (-2 * np.pi * b) - 1.0) < 1e-5
    ev = RA.valve_events(t, np.sin(2 * np.pi * t), 0.5)
    assert [o for _, o in ev] == [True, False]
    assert np.allclose([s for s, _ in ev], [1 / 12, 5 / 12], atol=1e-7)
