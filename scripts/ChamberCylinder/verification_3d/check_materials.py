# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Check that a generated 3D case and its 0D input describe the same material.

Evaluates the second Piola-Kirchhoff stress of both codes for random
deformation gradients and fiber directions, each as implemented:

  0D   ChamberCylinder compute_stress (src/model/ChamberCylinder.cpp) with
       the parameters of the 0D input: C1 (I1b - 3) + C2 (I2b - 3)
       + C3 exp(C4 (I1b - 3)^2) + C5 exp(C6 (I4b - 1)^2) plus the
       Holzapfel-Ogden terms (ho_*), as sums of dW/dIb dIb/dC, and the ST91
       pressure p = kappa/2 (J - 1/J)
  3D   svMultiPhysics compute_pk2cc (physics/solid/src/mat_models.cpp:
       neoHookean, MR, HGO, HolzapfelOgden as fictitious stress + bar_to_iso,
       ST91) with the parameters read from solver.xml

and compares the fiber and sheet fields written to the case (FIB_DIR of
fiber.vtu, sheet.vtu) with the 0D directions at their nodes or centroids.

Usage: python check_materials.py MODEL_0D.json CASE_DIR
"""
import argparse
import json
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np


def ho_family(E, a, b, k):
    """dW/dI of a Holzapfel-Ogden fiber or sheet term (ChamberCylinder.cpp)."""
    if a == 0.0:
        return 0.0
    o = 1.0 / (np.exp(k * E) + 1.0)
    rexp = np.exp(b * E * E)
    return a * ((1.0 - o) * E * rexp + 0.5 * k * (o - o * o) / b * (rexp - 1.0))


def stress_0d(F, f, s, v):
    """2nd PK stress of the ChamberCylinder material (incl. the ST91 pressure)."""
    C = F.T @ F
    J = np.linalg.det(F)
    Ci = np.linalg.inv(C)
    I = np.eye(3)
    M = np.outer(f, f)
    I1, I4 = np.trace(C), np.trace(C @ M)
    I2 = 0.5 * (I1 ** 2 - np.sum(C * C))
    Jm23 = J ** (-2.0 / 3.0)
    Jm43 = Jm23 * Jm23
    I1b, I4b = Jm23 * I1, Jm23 * I4
    w1 = v["C1"] + 2.0 * v["C3"] * v["C4"] * (I1b - 3.0) * np.exp(v["C4"] * (I1b - 3.0) ** 2)
    w2 = v["C2"]
    w4 = 2.0 * v["C5"] * v["C6"] * (I4b - 1.0) * np.exp(v["C6"] * (I4b - 1.0) ** 2)
    dI1b = Jm23 * (I - I1 / 3.0 * Ci)
    dI2b = Jm43 * (I1 * I - C - 2.0 / 3.0 * I2 * Ci)
    dI4b = Jm23 * (M - I4 / 3.0 * Ci)
    S = 2.0 * (w1 * dI1b + w2 * dI2b + w4 * dI4b)
    # Holzapfel-Ogden terms
    h = {k: v.get("ho_" + k, 0.0) for k in ("a", "b", "af", "bf", "as", "bs", "afs", "bfs", "k")}
    Hs, Hfs = np.outer(s, s), 0.5 * (np.outer(f, s) + np.outer(s, f))
    I6, I8 = np.sum(C * Hs), np.sum(C * Hfs)
    I6b, I8b = Jm23 * I6, Jm23 * I8
    dI6b, dI8b = Jm23 * (Hs - I6 / 3.0 * Ci), Jm23 * (Hfs - I8 / 3.0 * Ci)
    wa = 0.5 * h["a"] * np.exp(h["b"] * (I1b - 3.0))
    wf = ho_family(I4b - 1.0, h["af"], h["bf"], h["k"])
    ws = ho_family(I6b - 1.0, h["as"], h["bs"], h["k"])
    wfs = h["afs"] * I8b * np.exp(h["bfs"] * I8b ** 2)
    S += 2.0 * (wa * dI1b + wf * dI4b + ws * dI6b + wfs * dI8b)
    kappa = v.get("bulk_modulus", 0.0)
    return S + 0.5 * kappa * (J - 1.0 / J) * J * Ci


def svmp_material(xml_path):
    """Material parameters of an svMultiPhysics solver.xml (struct/ustruct); with
    several domains, those of domain 1 (the wall)."""
    eq = ET.parse(xml_path).getroot().find("Add_equation")
    val = lambda e, tag, d=None: float(e.findtext(tag)) if e.findtext(tag) else d  # noqa: E731
    dom = next((d for d in eq.findall("Domain") if d.get("id", "").strip() == "1"), eq)
    cm = dom.find("Constitutive_model")
    m = dict(equation=eq.get("type"), model=cm.get("type"), E=val(dom, "Elasticity_modulus"),
             nu=val(dom, "Poisson_ratio"), Kpen=val(dom, "Penalty_parameter", 0.0),
             vol=(dom.findtext("Dilational_penalty_model") or "").strip())
    for tag in ("c1", "c2", "a4", "b4", "a6", "b6", "kappa",
                "a", "b", "a4f", "b4f", "a4s", "b4s", "afs", "bfs", "k"):
        m[tag] = val(cm, tag, 0.0)
    return m


def stress_3d(F, f, s, m):
    """2nd PK stress of svMultiPhysics compute_pk2cc for struct with ST91."""
    C = F.T @ F
    J = np.linalg.det(F)
    Ci = np.linalg.inv(C)
    I = np.eye(3)
    J2d = J ** (-2.0 / 3.0)
    Inv1 = J2d * np.trace(C)
    mu = 0.5 * m["E"] / (1.0 + m["nu"])
    if m["model"] in ("neoHookean", "nHK"):
        Sb = 2.0 * (0.5 * mu) * I
    elif m["model"] in ("MR", "Mooney-Rivlin"):
        Sb = 2.0 * (m["c1"] + Inv1 * m["c2"]) * I - 2.0 * m["c2"] * J2d * C
    elif m["model"] == "HGO":
        k = m["kappa"]
        Hff, Hss = np.outer(f, f), np.outer(s, s)
        Inv4, Inv6 = J2d * f @ C @ f, J2d * s @ C @ s
        Eff = k * Inv1 + (1.0 - 3.0 * k) * Inv4 - 1.0
        Ess = k * Inv1 + (1.0 - 3.0 * k) * Inv6 - 1.0
        g2 = m["a4"] * Eff * np.exp(m["b4"] * Eff ** 2)
        g3 = m["a6"] * Ess * np.exp(m["b6"] * Ess ** 2)
        Sb = 2.0 * (0.5 * mu * I + g2 * (k * I + (1 - 3 * k) * Hff)
                    + g3 * (k * I + (1 - 3 * k) * Hss))
    elif m["model"] == "HolzapfelOgden":  # mat_models.cpp stIso_HO
        Hff, Hss = np.outer(f, f), np.outer(s, s)
        Hfs = 0.5 * (np.outer(f, s) + np.outer(s, f))
        Eff, Ess, Efs = J2d * f @ C @ f - 1.0, J2d * s @ C @ s - 1.0, J2d * f @ C @ s
        k = m["k"]
        of, os_ = 1.0 / (np.exp(k * Eff) + 1.0), 1.0 / (np.exp(k * Ess) + 1.0)
        c4f, c4s = 1.0 - of, 1.0 - os_
        dc4f, dc4s = k * (of - of ** 2), k * (os_ - os_ ** 2)
        Sb = m["a"] * np.exp(m["b"] * (Inv1 - 3.0)) * I
        Sb = Sb + 2.0 * m["afs"] * Efs * np.exp(m["bfs"] * Efs ** 2) * Hfs
        rexp = np.exp(m["b4f"] * Eff ** 2)
        Sb = Sb + 2.0 * m["a4f"] * (c4f * Eff * rexp + 0.5 * dc4f / m["b4f"] * (rexp - 1.0)) * Hff
        rexp = np.exp(m["b4s"] * Ess ** 2)
        Sb = Sb + 2.0 * m["a4s"] * (c4s * Ess * rexp + 0.5 * dc4s / m["b4s"] * (rexp - 1.0)) * Hss
    else:
        raise ValueError(f"unsupported model {m['model']}")
    r1 = J2d * np.sum(C * Sb) / 3.0
    S = J2d * Sb - r1 * Ci
    if m["vol"] == "ST91" and m["Kpen"] > 0.0:
        S += 0.5 * m["Kpen"] * (J - 1.0 / J) * J * Ci
    return S


def read_directions(path):
    """Locations (nodes or element centroids) and FIB_DIR of a fiber file."""
    from vtkmodules.util.numpy_support import vtk_to_numpy
    from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
    r = vtkXMLUnstructuredGridReader()
    r.SetFileName(path)
    r.Update()
    g = r.GetOutput()
    pts = vtk_to_numpy(g.GetPoints().GetData())
    conn = vtk_to_numpy(g.GetCells().GetConnectivityArray()).reshape(g.GetNumberOfCells(), -1)
    if g.GetPointData().GetArray("FIB_DIR") is not None:
        return pts, vtk_to_numpy(g.GetPointData().GetArray("FIB_DIR"))
    return pts[conn[:, :4]].mean(axis=1), vtk_to_numpy(g.GetCellData().GetArray("FIB_DIR"))


def check(model_0d, case_dir, n=2000, seed=0):
    """Max relative stress difference and max fiber-direction difference."""
    with open(model_0d) as fh:
        cfg = json.load(fh)
    v = next(b for b in cfg["vessels"] if b["zero_d_element_type"] == "ChamberCylinder")
    v = {k: float(x) for k, x in v["zero_d_element_values"].items()}
    for k in ("C1", "C2", "C3", "C4", "C5", "C6"):
        v.setdefault(k, 0.0)
    m = svmp_material(os.path.join(case_dir, "solver.xml"))
    if m["equation"] != "struct":
        raise ValueError("the stress check covers struct (penalty) cases")
    rng = np.random.default_rng(seed)
    err = 0.0
    for _ in range(n):
        F = np.eye(3) + rng.uniform(-0.3, 0.3, (3, 3))
        if np.linalg.det(F) < 0.2:
            continue
        f = rng.normal(size=3)
        f /= np.linalg.norm(f)
        s = np.cross(f, rng.normal(size=3))
        s /= np.linalg.norm(s)
        S0, S3 = stress_0d(F, f, s, v), stress_3d(F, f, s, m)
        err = max(err, np.abs(S0 - S3).max() / np.abs(S0).max())
    fib_err = None
    if os.path.exists(os.path.join(case_dir, "mesh", "fiber.vtu")):
        fib_err = 0.0
        for name in ("fiber", "sheet"):
            c, d = read_directions(os.path.join(case_dir, "mesh", f"{name}.vtu"))
            R, th = np.hypot(c[:, 0], c[:, 1]), np.arctan2(c[:, 1], c[:, 0])
            if name == "fiber":
                t = (R - v["Ri"]) / (v["Re"] - v["Ri"])
                a = np.radians(v["alpha_endo"] * (1 - t) + v["alpha_epi"] * t)
                d0 = np.stack([-np.sin(th) * np.cos(a), np.cos(th) * np.cos(a), np.sin(a)], axis=1)
            else:  # sheet e_R
                d0 = np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], axis=1)
            fib_err = max(fib_err, float(np.abs(d - d0).max()))
    return float(err), fib_err


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("model_0d")
    ap.add_argument("case_dir")
    a = ap.parse_args()
    err, fib = check(a.model_0d, a.case_dir)
    print(f"max relative difference of the 2nd PK stress, 0D vs 3D: {err:.2e}")
    if fib is not None:
        print(f"max difference of the fiber and sheet directions: {fib:.2e}")
    sys.exit(0 if err < 1e-12 and (fib is None or fib < 1e-12) else 1)


if __name__ == "__main__":
    main()
