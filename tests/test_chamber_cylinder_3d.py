"""Generator of the 3D svMultiPhysics counterpart of a ChamberCylinder model
(scripts/ChamberCylinder/verification_3d/cylinder_3d.py): mesh geometry, the
0D-to-3D mapping and the lumen-volume integral. svMultiPhysics is not run."""
import json
import math
import os
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
    names = {"endo", "epi", "base", "top"} | ({"x0", "y0"} if sector == 90.0 else set())
    assert set(m["faces"]) == names
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
    assert C3.generate(str(model), str(tmp_path / "a"), nr=1)["sector"] == 360.0
    with pytest.raises(ValueError):
        C3.generate(str(model), str(tmp_path / "b"), nr=1, sector=90.0)


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
