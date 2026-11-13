# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""3D svMultiPhysics counterpart of a ChamberCylinder 0D model.

Reads an svZeroDSolver input file that contains one ChamberCylinder block and
writes a matching svMultiPhysics case: a structured tetrahedral mesh of the
same thick-walled cylinder, its boundary faces, the solver input file and the
cavity-pressure load. Also runs the case (Docker image simvascular/solver or a
local executable) and computes the lumen volume from the results.

Model correspondence (0D -> 3D; the 0D "material", "active_stress" and
"formulation" objects select the components, see ChamberCylinder.h):
  geometry      Ri, Re, length                -> mesh
  material      genet23: C1 (neo-Hookean),    -> neoHookean, MR, or HGO with
                C2 (Mooney-Rivlin), C5, C6       a4 = 2 C5 C6, b4 = C6 (one
                fiber term; C3 = 0 required      fiber family, a6 = 0)
                holzapfel_ogden: a ... k      -> HolzapfelOgden a ... k
  fibers        alpha_endo, alpha_epi         -> FIB_DIR files (cylinder_fibers.py),
                sheet e_R                        fiber and sheet
  ends          fixed_length                  -> u_z = 0 on both end faces,
                                                 free in-plane sliding
                lids (pressure-loaded)        -> end plate (closed_end)
  load          inlet PRESSURE ramp           -> follower pressure on "endo"
  volumetric    st91: bulk_modulus kappa      -> ST91 penalty kappa (struct)
                incompressible                -> ustruct with nu = 0.5
  active        uniform_unsteady              -> UniformUnsteady (coupled_3d.py)

Mesh: n_r elements through the wall (the refinement parameter); the
circumferential and axial divisions follow from the element size
h = (Re - Ri) / n_r. A 90-degree sector uses symmetry planes (u_y = 0 on
theta = 0, u_x = 0 on theta = 90 deg); helical fibers break this symmetry and
need the full cylinder (and the full length). Its rigid-body modes (x, y
translation, rotation about z) are held by in-plane dashpots on both ends
(Robin, no stiffness): within a load step they act as a spring of stiffness
damping / dt, at rest they carry no force.
Quadratic tetrahedra (order 2) have mid-edge nodes on the exact cylindrical
surfaces.

Usage:
  python cylinder_3d.py generate MODEL_0D.json CASE_DIR --nr 4 [options]
  python cylinder_3d.py run CASE_DIR [--np 4] [--exe PATH]
  python cylinder_3d.py post CASE_DIR
"""
import argparse
import glob
import itertools
import json
import math
import os
import queue
import re
import shutil
import subprocess
import threading
import time

import numpy as np

MMHG = 133.322387415  # Pa; solver inputs are SI, reports use mmHg

# --------------------------------------------------------------------------- #
#  0D input
# --------------------------------------------------------------------------- #


# Formulation of the 0D block when its input omits a component (genet23);
# passive and active are the types of its material and active_stress objects.
FORMULATION = dict(passive="genet23", volumetric="incompressible", active="genet23", ends="lids")

# 0D material constants under the flat names used by these scripts
HO_FLAT = {"a": "ho_a", "b": "ho_b", "a4f": "ho_af", "b4f": "ho_bf", "a4s": "ho_as",
           "b4s": "ho_bs", "afs": "ho_afs", "bfs": "ho_bfs", "k": "ho_k"}
GENET23_FLAT = {"genet23": {k: k for k in ("C1", "C2", "C3", "C4", "C5", "C6")},
                "mooney_rivlin": {"W1": "C1", "W2": "C2"},
                "exponential": {"C0": "C3", "C1": "C4", "C2": "C5", "C3": "C6"}}


def formulation(vessel):
    """Formulation of a ChamberCylinder vessel of a 0D input: volumetric law
    and ends, and the types of its material (passive) and active stress
    (active)."""
    form = dict(FORMULATION)
    form.update({k: x for k, x in vessel.get("formulation", {}).items() if not k.startswith("_")})
    form["passive"] = "holzapfel_ogden" if vessel["material"]["type"] == "holzapfel_ogden" \
        else "genet23"
    form["active"] = vessel["active_stress"]["type"]
    return form


def flat_values(vessel):
    """zero_d_element_values of a ChamberCylinder vessel plus its material
    constants under flat names: C1-C6 (genet23 energy; mooney_rivlin and
    exponential map onto it) or ho_a ... ho_k (holzapfel_ogden)."""
    vals = dict(vessel["zero_d_element_values"])
    mat = vessel["material"]
    names = HO_FLAT if mat["type"] == "holzapfel_ogden" else GENET23_FLAT[mat["type"]]
    vals.update({flat: float(mat.get(key, 0.0)) for key, flat in names.items()})
    return vals


def read_0d(path):
    """ChamberCylinder values, formulation and the inlet pressure history of a
    0D input."""
    with open(path) as f:
        cfg = json.load(f)
    blocks = [v for v in cfg["vessels"] if v["zero_d_element_type"] == "ChamberCylinder"]
    if len(blocks) != 1:
        raise ValueError(f"{path}: expected one ChamberCylinder vessel, found {len(blocks)}")
    blk = blocks[0]
    inlet = blk["boundary_conditions"].get("inlet")
    bc = next((b for b in cfg["boundary_conditions"] if b["bc_name"] == inlet), None)
    if bc is None or bc["bc_type"] != "PRESSURE":
        raise ValueError(f"{path}: the ChamberCylinder inlet must be a PRESSURE boundary condition")
    return dict(name=blk["vessel_name"], values=flat_values(blk),
                formulation=formulation(blk), t=list(bc["bc_values"]["t"]),
                P=list(bc["bc_values"]["P"]), config=cfg)


def pressure_ramp(t, P):
    """(T_ramp, P_max) of a load P = 0 at t = 0, linear to P_max, then held."""
    i = int(np.argmax(P))
    t, P = np.asarray(t, float), np.asarray(P, float)
    if t[0] != 0.0 or P[0] != 0.0 or i == 0 or np.any(P[i:] != P[i]):
        raise ValueError("pressure must start at 0 at t = 0, ramp to its maximum and stay there")
    T, Pmax = t[i], P[i]
    if not np.allclose(P[: i + 1], Pmax * t[: i + 1] / T, rtol=1e-12, atol=0.0):
        raise ValueError("pressure must increase linearly up to its maximum")
    return float(T), float(Pmax)


HO_KEYS = ("a", "b", "af", "bf", "as", "bs", "afs", "bfs", "k")  # 0D ho_<key>
HO_SVMP = ("a", "b", "a4f", "b4f", "a4s", "b4s", "afs", "bfs", "k")  # svMultiPhysics


def material(values, form=None):
    """svMultiPhysics constitutive model equivalent to the 0D passive and
    volumetric laws of the formulation form (default: genet23)."""
    form = form or FORMULATION
    v = {k: float(values.get(k, 0.0))
         for k in ("C1", "C2", "C3", "C5", "C6", "gamma", "bulk_modulus")}
    ho = {k: float(values.get("ho_" + k, 0.0)) for k in HO_KEYS}
    if v["gamma"] != 0.0:
        raise ValueError("viscosity gamma has no counterpart in this 3D setup; set gamma = 0")
    if form["passive"] == "genet23" and v["C3"] != 0.0:
        raise ValueError("C3 exp(C4 (I1 - 3)^2) has no svMultiPhysics counterpart; set C3 = 0")
    # ST91 bulk modulus; 0 = incompressible
    kappa = v["bulk_modulus"] if form["volumetric"] == "st91" else 0.0
    fibers = dict(alpha_endo=float(values["alpha_endo"]), alpha_epi=float(values["alpha_epi"]))
    if form["passive"] == "holzapfel_ogden":
        if any(v[k] != 0.0 for k in ("C1", "C2", "C3", "C5")):
            raise ValueError("the passive law holzapfel_ogden takes no C1-C6")
        if ho["a"] <= 0.0:
            raise ValueError("HolzapfelOgden needs ho_a > 0 (small-strain shear modulus)")
        # mu = a is the small-strain shear modulus of the isotropic term.
        return dict(model="HolzapfelOgden", mu=ho["a"], kappa=kappa, **ho, **fibers)
    if any(x != 0.0 for x in ho.values()):
        raise ValueError('ho_* need the 0D material "holzapfel_ogden"')
    if v["C5"] != 0.0:
        # C5 exp(C6 (I4b - 1)^2) has the stress of the HGO fiber term
        # a4/(2 b4) (exp(b4 (I4b - 1)^2) - 1) with a4 = 2 C5 C6, b4 = C6.
        if v["C2"] != 0.0:
            raise ValueError("HGO has a neo-Hookean ground matrix; set C2 = 0 with fibers")
        return dict(model="HGO", mu=2.0 * v["C1"], a4=2.0 * v["C5"] * v["C6"], b4=v["C6"],
                    kappa=kappa, **fibers)
    if v["C2"] == 0.0:  # W = C1 (I1b - 3): svMultiPhysics neo-Hookean, C10 = mu / 2
        return dict(model="neoHookean", mu=2.0 * v["C1"], kappa=kappa)
    return dict(model="MR", mu=2.0 * (v["C1"] + v["C2"]), c1=v["C1"], c2=v["C2"], kappa=kappa)


# --------------------------------------------------------------------------- #
#  Mesh
# --------------------------------------------------------------------------- #

def helical(mat):
    """True for helical fibers: no mirror symmetry about theta = const or
    Z = const planes (circumferential or axial fibers have it)."""
    fibers = mat["model"] in ("HGO", "HolzapfelOgden")
    return fibers and not (mat["alpha_endo"] == mat["alpha_epi"] and mat["alpha_endo"] % 90.0 == 0.0)


# Six tetrahedra per hexahedron along its main diagonal (conforming on a
# structured grid). Corners are (di, dj, dk) offsets.
_KUHN = [[(0, 0, 0), tuple(int(a == p[0]) for a in range(3)),
          tuple(int(a in (p[0], p[1])) for a in range(3)), (1, 1, 1)]
         for p in itertools.permutations(range(3))]
_TET_EDGES = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]  # VTK quadratic tetra
_TRI_EDGES = [(0, 1), (1, 2), (2, 0)]                          # VTK quadratic triangle
_TET_FACES = [(0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)]


def mesh_cylinder(Ri, Re, L, nr, sector=90.0, order=2, ntheta=None, nz=None):
    """Structured tetrahedral mesh of the wall Ri <= R <= Re, 0 <= Z <= L.

    Returns points (N, 3), cells (M, 4 or 10), the coarse divisions and the
    boundary faces {name: (connectivity (K, 3 or 6), parent cell (K,))}, with
    face nodes ordered so the normal points out of the solid.
    """
    if sector not in (90.0, 360.0):
        raise ValueError("sector must be 90 (symmetry planes) or 360 degrees")
    h = (Re - Ri) / nr
    th = math.radians(sector)
    periodic = sector == 360.0
    nt = ntheta or max(1, math.ceil(th * 0.5 * (Ri + Re) / h))
    if periodic and not ntheta:  # node planes at theta = 0, 90, 180, 270 deg
        nt = 4 * math.ceil(nt / 4)
    nz = nz or max(1, math.ceil(L / h))
    s = 2 if order == 2 else 1

    index, points = {}, []

    def node(I, J, K):  # fine-grid indices; s * coarse index at the corners
        if periodic:
            J %= s * nt
        key = (I, J, K)
        if key not in index:
            R = Ri + (Re - Ri) * I / (s * nr)
            t = th * J / (s * nt)
            index[key] = len(points)
            points.append((R * math.cos(t), R * math.sin(t), L * K / (s * nz)))
        return index[key]

    cells, corners = [], []
    for i, j, k in itertools.product(range(nr), range(nt), range(nz)):
        for tet in _KUHN:
            c = [(i + a, j + b, k + d) for a, b, d in tet]
            x = [np.array(points[node(s * ci, s * cj, s * ck)]) for ci, cj, ck in c]
            if np.dot(np.cross(x[1] - x[0], x[2] - x[0]), x[3] - x[0]) < 0.0:
                c[1], c[2] = c[2], c[1]
            conn = [node(s * ci, s * cj, s * ck) for ci, cj, ck in c]
            if order == 2:
                conn += [node(*(c[a][m] + c[b][m] for m in range(3))) for a, b in _TET_EDGES]
            cells.append(conn)
            corners.append(c)

    pts = np.array(points)
    count = {}
    for e, c in enumerate(corners):
        for f in _TET_FACES:
            key = tuple(sorted((c[a][0], c[a][1] % nt if periodic else c[a][1], c[a][2])
                               for a in f))
            count.setdefault(key, []).append((e, f))

    faces = {}
    planes = [("endo", 0, 0), ("epi", 0, nr), ("base", 2, 0), ("top", 2, nz)]
    if not periodic:
        planes += [("y0", 1, 0), ("x0", 1, nt)]
    for key, owners in count.items():
        if len(owners) != 1:
            continue
        e, f = owners[0]
        name = next((n for n, ax, val in planes if all(v[ax] == val for v in key)), None)
        if name is None:
            raise RuntimeError(f"unclassified boundary face {key}")
        c = corners[e]
        loc = list(f)
        opp = ({0, 1, 2, 3} - set(f)).pop()
        x = [pts[cells[e][a]] for a in loc]
        if np.dot(np.cross(x[1] - x[0], x[2] - x[0]), pts[cells[e][opp]] - x[0]) > 0.0:
            loc[1], loc[2] = loc[2], loc[1]
        conn = [cells[e][a] for a in loc]
        if order == 2:
            conn += [cells[e][4 + _TET_EDGES.index(tuple(sorted((loc[a], loc[b]))))]
                     for a, b in _TRI_EDGES]
        faces.setdefault(name, ([], []))
        faces[name][0].append(conn)
        faces[name][1].append(e)
    faces = {n: (np.array(c), np.array(p)) for n, (c, p) in faces.items()}
    return dict(points=pts, cells=np.array(cells), faces=faces, order=order,
                nr=nr, ntheta=nt, nz=nz, sector=sector)


def mesh_closed_cylinder(Ri, Re, Lm, nr, plate, order=2, ntheta=None, nz=None, nz_plate=None,
                         n_in=None):
    """Structured tetrahedral mesh of a 90-degree sector of a tube closed at one
    end: the wall Ri <= R <= Re, 0 <= Z <= Lm (Z = 0 is the mid-length plane)
    and an end plate 0 <= R <= Re, Lm <= Z <= Lm + plate.

    The grid of mesh_cylinder is continued through the plate and inward to the
    axis (n_in elements over 0 <= R < Ri); the hexahedra at the axis collapse
    to wedges, of which the three non-degenerate tetrahedra remain. Cells of
    the wall have domain 1, cells of the plate domain 2. Boundary faces: endo
    (R = Ri), cap (the plate's inner face Z = Lm, R < Ri; both carry the
    cavity pressure), epi (R = Re), base (Z = 0), top (the plate's outer
    face) and the symmetry planes y0, x0.
    """
    h = (Re - Ri) / nr
    th = 0.5 * math.pi
    nt = ntheta or max(1, math.ceil(th * 0.5 * (Ri + Re) / h))
    nz = nz or max(1, math.ceil(Lm / h))
    nzp = nz_plate or max(1, math.ceil(plate / h))
    n_in = n_in or max(1, math.ceil(Ri / h))
    s = 2 if order == 2 else 1
    I_ax = -s * n_in  # fine radial index of the axis

    def radius(I):
        return Ri + (Re - Ri) * I / (s * nr) if I >= 0 else Ri * (1.0 + I / (s * n_in))

    def zcoord(K):
        return Lm * K / (s * nz) if K <= s * nz else Lm + plate * (K - s * nz) / (s * nzp)

    index, points = {}, []

    def node(I, J, K):  # fine-grid indices; all J coincide on the axis
        key = (I, 0 if I == I_ax else J, K)
        if key not in index:
            R, t = radius(I), th * key[1] / (s * nt)
            index[key] = len(points)
            points.append((R * math.cos(t), R * math.sin(t), zcoord(K)))
        return index[key]

    def mid(a, b):  # mid-edge node; an edge to the axis is radial at the other end's theta
        if a[0] == -n_in and b[0] != -n_in:
            a, b = b, a
        J = s * a[1] if b[0] == -n_in else a[1] + b[1]
        return node(a[0] + b[0], J, a[2] + b[2])

    cells, corner_ids, domains = [], [], []
    for i, j, k in itertools.product(range(-n_in, nr), range(nt), range(nz + nzp)):
        if i < 0 and k < nz:  # cavity
            continue
        for tet in _KUHN:
            c = [(i + a, j + b, k + d) for a, b, d in tet]
            if i == -n_in and any(c[p][0] == c[q][0] == -n_in and c[p][2] == c[q][2]
                                  for p, q in itertools.combinations(range(4), 2)):
                continue  # two corners on the same point of the axis
            ids = [node(s * ci, s * cj, s * ck) for ci, cj, ck in c]
            x = [np.array(points[n]) for n in ids]
            vol = np.dot(np.cross(x[1] - x[0], x[2] - x[0]), x[3] - x[0])
            if vol < 0.0:
                c[1], c[2], ids[1], ids[2] = c[2], c[1], ids[2], ids[1]
            if abs(vol) < 1e-12 * h ** 3:
                raise RuntimeError(f"degenerate tetrahedron in hexahedron {(i, j, k)}")
            conn = list(ids)
            if order == 2:
                conn += [mid(c[a], c[b]) for a, b in _TET_EDGES]
            cells.append(conn)
            corner_ids.append(ids)
            domains.append(1 if k < nz else 2)

    pts = np.array(points)
    cells = np.array(cells)
    count = {}
    for e, ids in enumerate(corner_ids):
        for f in _TET_FACES:
            count.setdefault(tuple(sorted(ids[a] for a in f)), []).append((e, f))

    tol = 1e-9 * Re
    Ztop = Lm + plate

    def classify(x):
        r = np.hypot(x[:, 0], x[:, 1])
        for name, ok in (("y0", np.abs(x[:, 1]) < tol), ("x0", np.abs(x[:, 0]) < tol),
                         ("base", np.abs(x[:, 2]) < tol), ("top", np.abs(x[:, 2] - Ztop) < tol),
                         ("epi", np.abs(r - Re) < tol),
                         ("endo", (np.abs(r - Ri) < tol) & (x[:, 2] < Lm + tol)),
                         ("cap", (np.abs(x[:, 2] - Lm) < tol) & (r < Ri + tol))):
            if np.all(ok):
                return name
        raise RuntimeError(f"unclassified boundary face {x.tolist()}")

    faces = {}
    for key, owners in count.items():
        if len(owners) != 1:
            continue
        e, f = owners[0]
        loc = list(f)
        opp = ({0, 1, 2, 3} - set(f)).pop()
        x = pts[cells[e][loc]]
        if np.dot(np.cross(x[1] - x[0], x[2] - x[0]), pts[cells[e][opp]] - x[0]) > 0.0:
            loc[1], loc[2] = loc[2], loc[1]
        conn = [cells[e][a] for a in loc]
        if order == 2:
            conn += [cells[e][4 + _TET_EDGES.index(tuple(sorted((loc[a], loc[b]))))]
                     for a, b in _TRI_EDGES]
        name = classify(pts[conn[:3]])
        faces.setdefault(name, ([], []))
        faces[name][0].append(conn)
        faces[name][1].append(e)
    faces = {n: (np.array(c), np.array(p)) for n, (c, p) in faces.items()}
    return dict(points=pts, cells=cells, faces=faces, order=order, nr=nr, ntheta=nt, nz=nz,
                sector=90.0, domains=np.array(domains), plate=plate, nz_plate=nzp, n_in=n_in)


# --------------------------------------------------------------------------- #
#  VTK XML (ASCII) writers
# --------------------------------------------------------------------------- #


def _array(name, data, vtype, ncomp=1):
    fmt = "%.17g" if vtype.startswith("Float") else "%d"
    rows = np.atleast_2d(np.asarray(data).reshape(len(data), -1))
    body = "\n".join(" ".join(fmt % v for v in r) for r in rows)
    nc = f' NumberOfComponents="{ncomp}"' if ncomp > 1 else ""
    return (f'<DataArray type="{vtype}" Name="{name}"{nc} format="ascii">\n'
            f"{body}\n</DataArray>\n")


def write_vtu(path, points, cells, vtk_type, point_data, cell_data):
    n = cells.shape[1]
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="1.0" '
                'byte_order="LittleEndian" header_type="UInt64">\n<UnstructuredGrid>\n'
                f'<Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(cells)}">\n')
        f.write("<PointData>\n" + "".join(_array(k, *v) for k, v in point_data.items())
                + "</PointData>\n")
        f.write("<CellData>\n" + "".join(_array(k, *v) for k, v in cell_data.items())
                + "</CellData>\n")
        f.write("<Points>\n" + _array("Points", points, "Float64", 3) + "</Points>\n")
        f.write("<Cells>\n" + _array("connectivity", cells.ravel(), "Int64")
                + _array("offsets", n * np.arange(1, len(cells) + 1), "Int64")
                + _array("types", np.full(len(cells), vtk_type), "UInt8") + "</Cells>\n")
        f.write("</Piece>\n</UnstructuredGrid>\n</VTKFile>\n")


def write_vtp(path, points, polys, point_data, cell_data):
    n = polys.shape[1]
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="PolyData" version="1.0" '
                'byte_order="LittleEndian" header_type="UInt64">\n<PolyData>\n'
                f'<Piece NumberOfPoints="{len(points)}" NumberOfVerts="0" NumberOfLines="0" '
                f'NumberOfStrips="0" NumberOfPolys="{len(polys)}">\n')
        f.write("<PointData>\n" + "".join(_array(k, *v) for k, v in point_data.items())
                + "</PointData>\n")
        f.write("<CellData>\n" + "".join(_array(k, *v) for k, v in cell_data.items())
                + "</CellData>\n")
        f.write("<Points>\n" + _array("Points", points, "Float64", 3) + "</Points>\n")
        f.write("<Polys>\n" + _array("connectivity", polys.ravel(), "Int64")
                + _array("offsets", n * np.arange(1, len(polys) + 1), "Int64") + "</Polys>\n")
        f.write("</Piece>\n</PolyData>\n</VTKFile>\n")


def face_extension(order):
    """Face file format: VTK PolyData has no quadratic triangle, so 6-node faces
    are written as unstructured grids (svMultiPhysics reads both)."""
    return "vtu" if order == 2 else "vtp"


def write_mesh(mesh, folder):
    """Volume mesh and faces in the SimVascular layout (1-based global IDs)."""
    os.makedirs(os.path.join(folder, "mesh-surfaces"), exist_ok=True)
    pts, cells = mesh["points"], mesh["cells"]
    quad = mesh["order"] == 2
    gid = np.arange(1, len(pts) + 1)
    cdata = {"GlobalElementID": (np.arange(1, len(cells) + 1), "Int32"),
             "ModelRegionID": (np.ones(len(cells), int), "Int32")}
    if "domains" in mesh:  # svMultiPhysics Domain_file_path reads DOMAIN_ID
        cdata["DOMAIN_ID"] = (mesh["domains"], "Int32")
    write_vtu(os.path.join(folder, "mesh-complete.mesh.vtu"), pts, cells, 24 if quad else 10,
              {"GlobalNodeID": (gid, "Int32")}, cdata)
    ext = face_extension(mesh["order"])
    for fid, (name, (conn, parent)) in enumerate(sorted(mesh["faces"].items()), start=1):
        used = np.unique(conn)
        local = np.searchsorted(used, conn)
        path = os.path.join(folder, "mesh-surfaces", f"{name}.{ext}")
        pdata = {"GlobalNodeID": (used + 1, "Int32")}
        cdata = {"GlobalElementID": (parent + 1, "Int32"),
                 "ModelFaceID": (np.full(len(conn), fid), "Int32")}
        if quad:
            write_vtu(path, pts[used], local, 22, pdata, cdata)
        else:
            write_vtp(path, pts[used], local, pdata, cdata)


# --------------------------------------------------------------------------- #
#  svMultiPhysics input
# --------------------------------------------------------------------------- #

FORMULATIONS = {
    # quadratic displacement, linear pressure, exactly incompressible
    "ustruct-th": dict(equation="ustruct", order=2),
    # linear velocity-pressure with VMS stabilization, exactly incompressible
    "ustruct-p1": dict(equation="ustruct", order=1),
    # quadratic displacement with an ST91 bulk penalty (nearly incompressible)
    "struct-p2": dict(equation="struct", order=2),
    # linear displacement with an ST91 bulk penalty (locks as kappa grows)
    "struct-p1": dict(equation="struct", order=1),
}


def solver_xml(case):
    m, f = case["material"], case["formulation"]
    eq = FORMULATIONS[f]["equation"]
    mu, kappa = m["mu"], m["kappa"]
    if kappa > 0.0:  # ST91 penalty as in the 0D block; nu, E only set C10 = mu / 2
        nu = (3.0 * kappa - 2.0 * mu) / (2.0 * (3.0 * kappa + mu))
        vol = (f"   <Dilational_penalty_model> ST91 </Dilational_penalty_model>\n"
               f"   <Penalty_parameter> {kappa:.17g} </Penalty_parameter>\n")
    elif eq == "struct":
        raise ValueError("struct needs a bulk modulus: set bulk_modulus > 0 in the 0D input")
    else:
        nu, vol = 0.5, ""
    E = 2.0 * mu * (1.0 + nu)
    if m["model"] == "MR":
        model = (f'   <Constitutive_model type="MR">\n      <c1> {m["c1"]:.17g} </c1>\n'
                 f'      <c2> {m["c2"]:.17g} </c2>\n   </Constitutive_model>\n')
    elif m["model"] == "HGO":  # second family (sheet file) switched off
        model = (f'   <Constitutive_model type="HGO">\n      <a4> {m["a4"]:.17g} </a4>\n'
                 f'      <b4> {m["b4"]:.17g} </b4>\n      <a6> 0.0 </a6>\n      <b6> 1.0 </b6>\n'
                 f'      <kappa> 0.0 </kappa>\n   </Constitutive_model>\n')
    elif m["model"] == "HolzapfelOgden":
        model = ('   <Constitutive_model type="HolzapfelOgden">\n'
                 + "".join(f"      <{t}> {m[k]:.17g} </{t}>\n" for k, t in zip(HO_KEYS, HO_SVMP))
                 + "   </Constitutive_model>\n")
    else:
        model = '   <Constitutive_model type="neoHookean"> </Constitutive_model>\n'
    extra = ""
    if f == "ustruct-th":
        extra = "   <Use_taylor_hood_type_basis> true </Use_taylor_hood_type_basis>\n"
    if f == "ustruct-p1":
        extra = ("   <Momentum_stabilization_coefficient> 1e-3 </Momentum_stabilization_coefficient>\n"
                 "   <Continuity_stabilization_coefficient> 1e-3 </Continuity_stabilization_coefficient>\n")
    ls = case["linear_solver"]
    assembly = (f"         <Assembly> {ls['assembly']} </Assembly>\n"
                if ls.get("assembly") else "")
    ext = face_extension(case["order"])
    faces = "".join(f'  <Add_face name="{n}">\n      <Face_file_path> mesh/mesh-surfaces/{n}.{ext} '
                    f"</Face_file_path>\n  </Add_face>\n" for n in case["faces"])
    tag = f"Fiber_direction_from_file_{case.get('fiber_location', 'cells')}"
    faces += "".join(f"  <{tag}> {fp} </{tag}>\n" for fp in case.get("fiber_files", []))
    # Only the nodal fields: svMultiPhysics' stress-type outputs (Stress,
    # Cauchy_stress, VonMises_stress, Strain, Def_grad, Jacobian,
    # Fiber_alignment, Volume) abort in get_nnx ("Error in computing shape
    # functions") for quadratic elements with mid-edge nodes on the curved
    # surfaces; post() computes them from the displacement instead.
    out = "     <Velocity> true </Velocity>\n     <Acceleration> true </Acceleration>\n"
    if case.get("fiber_files"):
        out += "     <Fiber_direction> true </Fiber_direction>\n"
    if eq == "ustruct":
        out += "     <Pressure> true </Pressure>\n"

    def dirichlet(name, direction):
        state = ("      <Impose_on_state_variable_integral> true "
                 "</Impose_on_state_variable_integral>\n" if eq == "ustruct" else "")
        return (f'   <Add_BC name="{name}" >\n      <Type> Dir </Type>\n      <Value> 0.0 </Value>\n'
                f"      <Effective_direction> {direction} </Effective_direction>\n{state}"
                f"   </Add_BC>\n")

    active = ""
    if case.get("active_tension_file"):  # svMultiPhysics UniformUnsteady
        active = ("   <Active_stress>\n      <Model> UniformUnsteady </Model>\n"
                  "      <Directional_distribution>\n"
                  "         <Fiber_direction> 1.0 </Fiber_direction>\n"
                  "         <Sheet_direction> 0.0 </Sheet_direction>\n"
                  "         <Sheet_normal_direction> 0.0 </Sheet_normal_direction>\n"
                  "      </Directional_distribution>\n      <UniformUnsteady>\n"
                  "         <Ramp> false </Ramp>\n"
                  f"         <Temporal_values_file_path> {case['active_tension_file']} "
                  "</Temporal_values_file_path>\n      </UniformUnsteady>\n   </Active_stress>\n")
    plate = case.get("plate", 0.0) > 0.0
    material = (f"   <Density> {case['density']:.17g} </Density>\n"
                f"   <Elasticity_modulus> {E:.17g} </Elasticity_modulus>\n"
                f"   <Poisson_ratio> {nu} </Poisson_ratio>\n{model}{vol}{active}")
    domain_file = ""
    if plate:  # wall (domain 1) and an isotropic end plate (domain 2)
        if active or case.get("coupling"):
            raise ValueError("the closed end is implemented for passive inflation only")
        pm = case["plate_material"]
        nu_p = (3.0 * pm["kappa"] - 2.0 * pm["mu"]) / (2.0 * (3.0 * pm["kappa"] + pm["mu"]))
        indent = lambda s: "".join("   " + ln + "\n" for ln in s.splitlines())  # noqa: E731
        material = ('   <Domain id="1" >\n' + indent(material) + "   </Domain>\n"
                    '   <Domain id="2" >\n'
                    f"      <Density> {case['density']:.17g} </Density>\n"
                    f"      <Elasticity_modulus> {2.0 * pm['mu'] * (1.0 + nu_p):.17g} </Elasticity_modulus>\n"
                    f"      <Poisson_ratio> {nu_p} </Poisson_ratio>\n"
                    '      <Constitutive_model type="neoHookean"> </Constitutive_model>\n'
                    "      <Dilational_penalty_model> ST91 </Dilational_penalty_model>\n"
                    f"      <Penalty_parameter> {pm['kappa']:.17g} </Penalty_parameter>\n"
                    "   </Domain>\n")
        domain_file = "  <Domain_file_path> mesh/mesh-complete.mesh.vtu </Domain_file_path>\n"
    cp = case.get("coupling")
    if cp:  # cavity pressure from svZeroDSolver, flow of the endo face to it
        endo = ("   <Add_BC name=\"endo\" >\n      <Type> Neu </Type>\n"
                "      <Time_dependence> Coupled </Time_dependence>\n"
                "      <Follower_pressure_load> true </Follower_pressure_load>\n"
                "      <Coupling_interface>\n"
                f"         <svZeroDSolver_block> {cp['block']} </svZeroDSolver_block>\n"
                "      </Coupling_interface>\n   </Add_BC>\n")
        interface = ("   <svZeroDSolver_interface>\n"
                     "      <Coupling_type> implicit </Coupling_type>\n"
                     f"      <Configuration_file> {cp['config']} </Configuration_file>\n"
                     f"      <Shared_library> {cp['library']} </Shared_library>\n"
                     "      <Initial_flows> 0.0 </Initial_flows>\n"
                     f"      <Initial_pressures> {cp['initial_pressure']:.17g} </Initial_pressures>\n"
                     # finite-difference step of dP/dQ: max(rms(Q) * rel, abs);
                     # the defaults (1e-5, 1e-7) are ~10% of SI cavity flows
                     "      <Finite_difference_absolute_perturbation> 1e-14 "
                     "</Finite_difference_absolute_perturbation>\n"
                     "      <Finite_difference_relative_perturbation> 1e-8 "
                     "</Finite_difference_relative_perturbation>\n"
                     "   </svZeroDSolver_interface>\n")
    else:  # the closed end's plate carries the cavity pressure on its inner face
        endo = "".join(f"   <Add_BC name=\"{n}\" >\n      <Type> Neu </Type>\n"
                       "      <Time_dependence> Unsteady </Time_dependence>\n"
                       "      <Temporal_values_file_path> endo_pressure.dat </Temporal_values_file_path>\n"
                       "      <Ramp_function> true </Ramp_function>\n"
                       "      <Follower_pressure_load> true </Follower_pressure_load>\n   </Add_BC>\n"
                       for n in (("endo", "cap") if plate else ("endo",)))
        interface = ""
    # u_z = 0 on the mid-length plane (base) and, with fixed length, on the end
    bcs = dirichlet("base", "(0, 0, 1)") + ("" if plate else dirichlet("top", "(0, 0, 1)"))
    if "y0" in case["faces"]:  # symmetry planes of the sector
        bcs += dirichlet("y0", "(0, 1, 0)") + dirichlet("x0", "(1, 0, 0)")
    if case.get("end_damping", 0.0) > 0.0:  # full cylinder: in-plane dashpots on both ends
        bcs += "".join(f'   <Add_BC name="{n}" >\n      <Type> Robin </Type>\n'
                       "      <Stiffness> 0.0 </Stiffness>\n"
                       f"      <Damping> {case['end_damping']:.6g} </Damping>\n   </Add_BC>\n"
                       for n in case.get("dashpot_faces", ("base", "top")))
    if case.get("robin_stiffness", 0.0) > 0.0:
        bcs += (f'   <Add_BC name="epi" >\n      <Type> Robin </Type>\n'
                f"      <Stiffness> {case['robin_stiffness']:.6g} </Stiffness>\n"
                f"      <Damping> 0.0 </Damping>\n   </Add_BC>\n")
    return f"""<?xml version="1.0" encoding="UTF-8" ?>
<svMultiPhysicsFile version="0.1">
<!-- Generated by cylinder_3d.py from {case['source_0d']} (SI units) -->

<GeneralSimulationParameters>
  <Continue_previous_simulation> 0 </Continue_previous_simulation>
  <Number_of_spatial_dimensions> 3 </Number_of_spatial_dimensions>
  <Number_of_time_steps> {case['n_steps']} </Number_of_time_steps>
  <Time_step_size> {case['dt']:.17g} </Time_step_size>
  <Spectral_radius_of_infinite_time_step> 0.0 </Spectral_radius_of_infinite_time_step>
  <Searched_file_name_to_trigger_stop> STOP_SIM </Searched_file_name_to_trigger_stop>
  <Save_results_to_VTK_format> 1 </Save_results_to_VTK_format>
  <Name_prefix_of_saved_VTK_files> result </Name_prefix_of_saved_VTK_files>
  <Save_results_in_folder> results </Save_results_in_folder>
  <Increment_in_saving_VTK_files> 1 </Increment_in_saving_VTK_files>
  <Start_saving_after_time_step> 1 </Start_saving_after_time_step>
  <Increment_in_saving_restart_files> 1 </Increment_in_saving_restart_files>
  <Overwrite_restart_file> true </Overwrite_restart_file>
  <Convert_BIN_to_VTK_format> 0 </Convert_BIN_to_VTK_format>
  <Verbose> 1 </Verbose>
  <Warning> 1 </Warning>
  <Debug> 0 </Debug>
</GeneralSimulationParameters>

<Add_mesh name="msh" >
  <Mesh_file_path> mesh/mesh-complete.mesh.vtu </Mesh_file_path>
{faces}{domain_file}  <Mesh_scale_factor> 1.0 </Mesh_scale_factor>
</Add_mesh>

<Add_equation type="{eq}" >
   <Coupled> true </Coupled>
   <Min_iterations> 3 </Min_iterations>
   <Max_iterations> 20 </Max_iterations>
   <Tolerance> {case['newton_tol']:.3g} </Tolerance>
{extra}
{material}
   <Output type="Spatial" >
     <Displacement> true </Displacement>
{out}   </Output>

   <LS type="{ls['type']}" >
      <Linear_algebra type="{ls['algebra']}" >
         <Preconditioner> {ls['preconditioner']} </Preconditioner>
{assembly}      </Linear_algebra>
      <Tolerance> {ls['tol']:.3g} </Tolerance>
      <Absolute_tolerance> {ls['abs_tol']:.3g} </Absolute_tolerance>
      <Max_iterations> {ls['max_iter']} </Max_iterations>
      <Krylov_space_dimension> {ls['krylov']} </Krylov_space_dimension>
   </LS>

{interface}{endo}{bcs}
</Add_equation>
</svMultiPhysicsFile>
"""


def generate(model_0d, case_dir, nr, formulation="struct-p2", sector=None, ntheta=None,
             nz=None, n_steps=20, density=1.0, newton_tol=1e-10, linear_solver=None,
             robin_factor=0.0, half_length=False, fiber_location="points", end_damping=1e-3,
             dashpot_faces=("base", "top"), closed_end=False, plate_factor=1.0):
    """Write the svMultiPhysics case of a ChamberCylinder 0D input.

    half_length meshes 0 <= Z <= L/2 with u_z = 0 on Z = L/2, the mid-length
    symmetry plane of the Z-homogeneous solution (not with helical fibers:
    the plane would be a slip end). fiber_location "points" writes nodal
    fiber directions, "cells" one direction per element. The full cylinder
    gets in-plane dashpots on both ends against its rigid-body modes, with
    damping end_damping * mu / Re * dt: a spring of end_damping * mu / Re
    within a load step, no force at rest (dashpot_faces: the faces that carry
    them). robin_factor > 0 adds an epicardial spring of stiffness
    robin_factor * mu / Re.

    closed_end (0D formulation "ends": "lids", the pressure-loaded lids): the eighth
    model of the tube, 0 <= Z <= L/2 with the mid-length plane at Z = 0,
    closed at Z = L/2 by an isotropic plate of thickness plate_factor (Re -
    Ri) (neo-Hookean with the small-strain shear and bulk moduli of the
    wall) whose inner face carries the cavity pressure. The length is free:
    by equilibrium the wall then carries the axial force P pi r_i^2 of the 0D
    lids away from the plate, and the mid-length section is compared.
    """
    m0 = read_0d(model_0d)
    v, form = m0["values"], m0["formulation"]
    if closed_end and form["ends"] != "lids":
        raise ValueError('a closed end loads the wall axially: set the 0D formulation "ends": "lids"')
    if not closed_end and form["ends"] != "fixed_length":
        raise ValueError('the 3D model has fixed-length ends: set the 0D formulation '
                         '"ends": "fixed_length"')
    T, Pmax = pressure_ramp(m0["t"], m0["P"])
    mat = material(v, form)
    fibers = mat["model"] in ("HGO", "HolzapfelOgden")
    symmetric = not helical(mat)
    if sector is None:
        sector = 90.0 if symmetric else 360.0
    if not symmetric and sector != 360.0:
        raise ValueError("helical fibers break the mirror symmetry: use sector = 360")
    if not symmetric and half_length:
        raise ValueError("helical fibers: the mid-length plane would be a slip end, mesh the full length")
    if closed_end and not (symmetric and sector == 90.0):
        raise ValueError("the closed end is implemented for the eighth model (no helical fibers)")
    fm = FORMULATIONS[formulation]
    Lm = 0.5 * v["length"] if half_length or closed_end else v["length"]
    if closed_end:
        plate = plate_factor * (v["Re"] - v["Ri"])
        mesh = mesh_closed_cylinder(v["Ri"], v["Re"], Lm, nr, plate, fm["order"], ntheta, nz)
    else:
        mesh = mesh_cylinder(v["Ri"], v["Re"], Lm, nr, sector, fm["order"], ntheta, nz)
    write_mesh(mesh, os.path.join(case_dir, "mesh"))
    fiber_files = []
    if fibers:
        from cylinder_fibers import cylinder_fibers, write_fibers
        with np.errstate(invalid="ignore", divide="ignore"):  # the plate's axis nodes
            fib = cylinder_fibers(mesh["points"], mesh["cells"], v["Ri"], v["Re"],
                                  mat["alpha_endo"], mat["alpha_epi"], fiber_location)
        for k, e in (("fiber", 0), ("sheet", 1), ("normal", 2)):  # unused: isotropic plate
            bad = ~np.all(np.isfinite(fib[k]), axis=1)
            fib[k][bad] = np.eye(3)[e]
        paths = write_fibers(os.path.join(case_dir, "mesh"), mesh["points"], mesh["cells"],
                             24 if fm["order"] == 2 else 10, fib, fiber_location)
        fiber_files = [os.path.relpath(p, case_dir) for p in paths[:2]]  # fiber, sheet
    with open(os.path.join(case_dir, "endo_pressure.dat"), "w") as f:
        f.write(f"2 1\n0.0 0.0\n{T:.17g} {Pmax:.17g}\n")
    # BiCGStab with the iteration cap of the svMultiPhysics struct tests: the
    # solves stop early (inexact Newton) but Newton still converges in 4-5
    # iterations, 6x faster than GMRES to 1e-12. SI residuals are small
    # (~1e-5 N): the default absolute tolerance (1e-10) would stop the linear
    # solver before Newton converges.
    ls = dict(type="BICG", algebra="fsils", preconditioner="fsils", tol=1e-12,
              abs_tol=1e-20, max_iter=600, krylov=100)
    ls.update(linear_solver or {})
    case = dict(source_0d=os.path.basename(model_0d), block=m0["name"],
                Ri=v["Ri"], Re=v["Re"], length=v["length"], mesh_length=Lm, material=mat,
                fiber_files=fiber_files, fiber_location=fiber_location,
                robin_stiffness=robin_factor * mat["mu"] / v["Re"] if sector == 360.0 else 0.0,
                **({"end_damping": end_damping * mat["mu"] / v["Re"] * T / n_steps}
                   if sector == 360.0 else {}),
                **({"dashpot_faces": list(dashpot_faces)}
                   if sector == 360.0 and tuple(dashpot_faces) != ("base", "top") else {}),
                **(dict(plate=mesh["plate"], nz_plate=mesh["nz_plate"], n_in=mesh["n_in"],
                        plate_material=dict(mu=mat["mu"], kappa=mat["kappa"]))
                   if closed_end else {}),
                formulation=formulation, order=fm["order"], sector=sector, nr=nr,
                ntheta=mesh["ntheta"], nz=mesh["nz"], n_nodes=len(mesh["points"]),
                n_cells=len(mesh["cells"]), faces=sorted(mesh["faces"]),
                T_ramp=T, P_max=Pmax, n_steps=n_steps, dt=T / n_steps, density=density,
                newton_tol=newton_tol, linear_solver=ls)
    with open(os.path.join(case_dir, "solver.xml"), "w") as f:
        f.write(solver_xml(case))
    with open(os.path.join(case_dir, "case.json"), "w") as f:
        json.dump(case, f, indent=2)
    return case


# --------------------------------------------------------------------------- #
#  Run
# --------------------------------------------------------------------------- #

IMAGE = "simvascular/solver:latest"
EXE_IN_IMAGE = "/build-trilinos/svMultiPhysics-build/bin/svmultiphysics"


def docker_cli():
    """The docker executable (Docker Desktop on macOS installs it in ~/.docker/bin)."""
    for path in (shutil.which("docker"), os.path.expanduser("~/.docker/bin/docker"),
                 "/Applications/Docker.app/Contents/Resources/bin/docker"):
        if path and os.path.exists(path):
            return path
    raise FileNotFoundError("docker not found; install Docker or pass --exe")


class SolverFailure(RuntimeError):
    """svMultiPhysics crashed, diverged or exceeded its time limit."""


_NEWTON_ROW = re.compile(r"^\s*ST\s+(\d+)-(\d+)(s?)\s")


def run(case_dir, np_=1, exe=None, image=IMAGE, platform="linux/amd64", max_time=None,
        max_residual=1e6, echo=True, mounts=None, echo_every=1):
    """Run svMultiPhysics on a case: locally if exe is given, else in Docker.

    The solver output is streamed to CASE_DIR/svmultiphysics.log and checked
    line by line: the run is stopped and SolverFailure raised as soon as a
    Newton iteration reports NaN or a relative residual above max_residual, a
    time step ends unconverged (above 100x the Newton tolerance), the solver
    exits with an error, or max_time seconds pass. With echo, one line per
    echo_every time steps reports the Newton iterations and the final
    residual. mounts {host_dir: container_dir} adds Docker volumes.
    """
    case_dir = os.path.abspath(case_dir)
    with open(os.path.join(case_dir, "case.json")) as f:
        case = json.load(f)
    mpi = ["mpirun", "-n", str(np_)] if np_ > 1 else []
    env, name = None, None
    if exe:
        cmd = mpi + [exe, "solver.xml"]
    else:
        docker = docker_cli()
        env = dict(os.environ, PATH=os.path.dirname(docker) + os.pathsep + os.environ["PATH"])
        name = f"svmp-{os.getpid()}-{int(time.time() * 1e3)}"
        cmd = ([docker, "run", "--rm", "--name", name, "--platform", platform,
                "-e", "OMPI_ALLOW_RUN_AS_ROOT=1", "-e", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
                "-v", f"{case_dir}:/case", "-w", "/case"]
               + [a for h, c in (mounts or {}).items() for a in ("-v", f"{h}:{c}")] + [image]
               + (mpi[:1] + ["--allow-run-as-root"] + mpi[1:] if mpi else [])
               + [EXE_IN_IMAGE, "solver.xml"])
    tic = time.time()
    proc = subprocess.Popen(cmd, cwd=case_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, env=env, bufsize=1)
    lines = queue.Queue()
    threading.Thread(target=lambda: [lines.put(ln) for ln in proc.stdout] + [lines.put(None)],
                     daemon=True).start()
    check = LogChecker(case, max_residual, echo, echo_every)
    reason = None
    with open(os.path.join(case_dir, "svmultiphysics.log"), "w") as log:
        while reason is None:
            try:
                line = lines.get(timeout=1.0)
            except queue.Empty:
                line = ""
            if line is None:
                break
            log.write(line)
            log.flush()
            reason = check(line)
            if max_time and time.time() - tic > max_time:
                reason = f"time limit of {max_time:.0f} s reached"
    if reason is not None:  # stop the solver
        if name:
            subprocess.run([docker, "stop", "-t", "2", name], env=env, capture_output=True)
        proc.terminate()
    rc = proc.wait()
    if reason is None and rc != 0:
        reason = f"svMultiPhysics exited with code {rc}"
    if reason is not None:
        raise SolverFailure(f"{case_dir}: {reason}\n" + "\n".join(check.tail))
    return time.time() - tic


class LogChecker:
    """Line-by-line check of the svMultiPhysics output (see `run`): returns a
    failure reason for a Newton iteration with NaN or a relative residual above
    max_residual, or a time step that ends above 100x the Newton tolerance;
    with echo, prints one progress line per echo_every time steps. With
    strict=False a time step that ends above the tolerance is not a failure
    (svMultiPhysics goes on); it is listed in `unconverged` instead."""

    def __init__(self, case, max_residual=1e6, echo=True, echo_every=1, strict=True):
        self.case, self.max_residual = case, max_residual
        self.echo, self.echo_every, self.strict = echo, echo_every, strict
        self.tail, self.its, self.step, self.tic = [], 0, 0, time.time()
        self.unconverged = []  # (time step, final residual), strict=False

    def __call__(self, line):
        self.tail = (self.tail + [line.rstrip()])[-8:]
        m = _NEWTON_ROW.match(line)
        if not m:
            return None
        self.its += 1
        nums = re.findall(r"\d\.\d+e[+-]\d+|nan", line)
        res = float(nums[1]) if len(nums) > 1 else float("nan")  # Ri/R1 (step)
        res0 = float(nums[2]) if len(nums) > 2 else res  # Ri/R0 (first step)
        step = self.step = int(m.group(1))
        if not math.isfinite(res) or res > self.max_residual:
            return f"time step {step}, Newton iteration {m.group(2)}: residual {res:.1e}"
        if m.group(3):  # last iteration of the time step
            if min(res, res0) > 100.0 * self.case["newton_tol"]:
                if self.strict:
                    return f"time step {step}: Newton stopped at residual {res:.1e}"
                self.unconverged.append((step, min(res, res0)))
            n = self.case["n_steps"]
            if self.echo and (step % self.echo_every == 0 or step == n):
                print(f"    step {step:4d}/{n}: {self.its:2d} Newton its, residual {res:.1e}, "
                      f"{time.time() - self.tic:7.1f} s", flush=True)
            self.its = 0
        return None


# --------------------------------------------------------------------------- #
#  Post-processing
# --------------------------------------------------------------------------- #

def _tri_rule():
    """Degree-5 seven-point rule on the reference triangle: (xi, eta, weight)."""
    a1, w1 = 0.470142064105115, 0.132394152788506
    a2, w2 = 0.101286507323456, 0.125939180544827
    pts = [(1 / 3, 1 / 3, 0.225)]
    for a, w in ((a1, w1), (a2, w2)):
        pts += [(a, a, w), (1 - 2 * a, a, w), (a, 1 - 2 * a, w)]
    return [(x, y, 0.5 * w) for x, y, w in pts]  # weights sum to the area 1/2


def _tri6_shape(x, y):
    l = (1 - x - y, x, y)
    N = [l[0] * (2 * l[0] - 1), l[1] * (2 * l[1] - 1), l[2] * (2 * l[2] - 1),
         4 * l[0] * l[1], 4 * l[1] * l[2], 4 * l[2] * l[0]]
    dl = ((-1, -1), (1, 0), (0, 1))
    dN = [[(4 * l[i] - 1) * dl[i][d] for d in range(2)] for i in range(3)]
    for (i, j) in _TRI_EDGES:
        dN.append([4 * (dl[i][d] * l[j] + l[i] * dl[j][d]) for d in range(2)])
    return np.array(N), np.array(dN)


def lumen_volume(x, conn, sector):
    """Lumen volume enclosed by the deformed endo face between the end planes.

    V = int (x n_x + y n_y) / 2 dA over the endo face with n pointing into the
    wall (the end planes and the symmetry planes contribute nothing), scaled
    from the sector to the full cylinder.
    """
    V = 0.0
    if conn.shape[1] == 3:
        a, b, c = x[conn[:, 0]], x[conn[:, 1]], x[conn[:, 2]]
        n = 0.5 * np.cross(b - a, c - a)
        g = (a + b + c) / 3.0
        V = -0.5 * np.sum(g[:, 0] * n[:, 0] + g[:, 1] * n[:, 1])
    else:
        for xi, eta, w in _tri_rule():
            N, dN = _tri6_shape(xi, eta)
            xe = x[conn]                      # (K, 6, 3)
            p = np.einsum("a,kai->ki", N, xe)
            t1 = np.einsum("a,kai->ki", dN[:, 0], xe)
            t2 = np.einsum("a,kai->ki", dN[:, 1], xe)
            n = np.cross(t1, t2)
            V -= 0.5 * w * np.sum(p[:, 0] * n[:, 0] + p[:, 1] * n[:, 1])
    return V * 360.0 / sector


def _read_vtu(path):
    """Points and the Displacement field of an svMultiPhysics result file."""
    try:
        from vtkmodules.util.numpy_support import vtk_to_numpy
        from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
    except ImportError as e:  # pragma: no cover
        raise ImportError("post-processing needs VTK (pip install vtk)") from e
    r = vtkXMLUnstructuredGridReader()
    r.SetFileName(path)
    r.Update()
    g = r.GetOutput()
    pts = vtk_to_numpy(g.GetPoints().GetData())
    return pts, {"Displacement": vtk_to_numpy(g.GetPointData().GetArray("Displacement"))}


def kinematics(X, x, levels, z0=None, dz=None):
    """0D kinematic quantities of a deformed 3D cylinder.

    Fits the rotation of every node about the axis as beta Z + phi(R) and
    averages over each radial node level R_k (index `levels`) the radial
    displacement rho = r - R, the transmural rotation phi and the axial
    displacement u_z, the latter two relative to the endocardium: the 0D
    fields rho(R), phi(R) - phi(Ri) and eta(R) - eta(Ri). Also returns the
    largest deviation of rho from its level mean (non-homogeneity). With z0,
    only the nodes with |Z - z0| <= dz enter, and the rotation is fitted as
    phi(R) + beta (Z - z0) + gamma (Z - z0)^2: beta is the local twist rate
    at z0.
    """
    dth = np.arctan2(x[:, 1], x[:, 0]) - np.arctan2(X[:, 1], X[:, 0])
    dth = (dth + np.pi) % (2.0 * np.pi) - np.pi
    rho = np.hypot(x[:, 0], x[:, 1]) - np.hypot(X[:, 0], X[:, 1])
    uz = x[:, 2] - X[:, 2]
    nlev = levels.max() + 1
    Z = X[:, 2]
    if z0 is None:
        use, cols = np.ones(len(X), bool), [Z]
    else:
        use = np.abs(Z - z0) <= dz * (1.0 + 1e-9)
        cols = [Z - z0, (Z - z0) ** 2]
    lev = levels[use]
    if np.any(np.bincount(lev, minlength=nlev) == 0):
        raise ValueError("a radial node level has no nodes in the compared section")
    # rotation = beta Z (+ gamma Z^2) + phi_k: least squares, one offset per level
    nc = len(cols)
    A = np.zeros((use.sum(), nc + nlev))
    for j, c in enumerate(cols):
        A[:, j] = c[use]
    A[np.arange(use.sum()), nc + lev] = 1.0
    coef, *_ = np.linalg.lstsq(A, dth[use], rcond=None)
    mean = lambda f: np.bincount(lev, f[use]) / np.bincount(lev)  # noqa: E731
    rho_k = mean(rho)
    # axial strain: u_z = eps (Z - z0) + offset per level (eps = 0 for fixed length)
    B = np.zeros((use.sum(), 1 + nlev))
    B[:, 0] = Z[use] - (z0 or 0.0)
    B[np.arange(use.sum()), 1 + lev] = 1.0
    eps = np.linalg.lstsq(B, uz[use], rcond=None)[0][0]
    return dict(beta=float(coef[0]), rho=rho_k, phi=coef[nc:] - coef[nc],
                eta=mean(uz) - mean(uz)[0], eps=float(eps),
                rho_spread=float(np.max(np.abs(rho[use] - rho_k[lev]))))


def cylinder_fibers_at_nodes(X, Ri, Re, alpha_endo, alpha_epi):
    """Reference fiber and sheet (e_R) directions at the nodes (0D convention);
    unit Cartesian vectors on the axis (closed end), where they are undefined."""
    from cylinder_fibers import cylinder_fibers
    with np.errstate(invalid="ignore", divide="ignore"):
        fib = cylinder_fibers(X, None, Ri, Re, alpha_endo, alpha_epi, "points")
    for k, e in (("fiber", 0), ("sheet", 1), ("normal", 2)):
        fib[k][~np.all(np.isfinite(fib[k]), axis=1)] = np.eye(3)[e]
    return fib


def _tet_shape_derivatives(order):
    """dN_b/dxi (natural coordinates) of the VTK tetrahedron at its own nodes a,
    as an array (n, n, 3) indexed [a, b, k]."""
    I4 = np.eye(4)
    L_nodes = [I4[i] for i in range(4)]
    if order == 2:
        L_nodes += [0.5 * (I4[i] + I4[j]) for i, j in _TET_EDGES]
    dL = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    dN = np.zeros((len(L_nodes), len(L_nodes), 3))
    for a, L in enumerate(L_nodes):
        if order == 1:
            dN[a] = dL
            continue
        for i in range(4):
            dN[a, i] = (4.0 * L[i] - 1.0) * dL[i]
        for m, (i, j) in enumerate(_TET_EDGES):
            dN[a, 4 + m] = 4.0 * (L[i] * dL[j] + L[j] * dL[i])
    return dN


def nodal_deformation_gradient(X, cells, u, order):
    """F = I + du/dX at the nodes with the (curved) element geometry, averaged
    over the elements that share a node."""
    dN = _tet_shape_derivatives(order)
    Jx = np.einsum("abk,ebi->eaik", dN, X[cells])          # dX_i/dxi_k
    dNdX = np.einsum("abk,eaki->eabi", dN, np.linalg.inv(Jx))
    G = np.einsum("ebi,eabj->eaij", u[cells], dNdX)         # du_i/dX_j
    F = np.zeros((len(X), 3, 3))
    count = np.zeros(len(X))
    np.add.at(F, cells, G)
    np.add.at(count, cells, 1.0)
    return np.eye(3) + F / count[:, None, None]


_SYM = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0)]  # svMultiPhysics order


def _sym6(T):
    return np.stack([T[:, i, j] for i, j in _SYM], axis=1)


def struct_fields(X, cells, u, order, material, fiber, sheet):
    """The svMultiPhysics struct output fields at the nodes, from the
    displacement: deformation gradient, Jacobian, Green-Lagrange strain, 2nd
    Piola-Kirchhoff and Cauchy stress (svMultiPhysics material model of
    solver.xml, see check_materials.stress_3d), von Mises stress, fiber
    stretch and deformed fiber direction, plus the cylindrical components
    (rr, tt, zz, rt, tz, zr) of the strain (reference basis) and of the Cauchy
    stress (deformed basis)."""
    from check_materials import stress_3d
    out = kinematic_fields(X, cells, u, order, fiber)
    F = out["Def_grad"].reshape(-1, 3, 3)
    J = out["Jacobian"]
    S = np.array([stress_3d(F[n], fiber[n], sheet[n], material) for n in range(len(X))])
    sig = np.einsum("nik,nkl,njl->nij", F, S, F) / J[:, None, None]
    dev = sig - np.trace(sig, axis1=1, axis2=2)[:, None, None] * np.eye(3) / 3.0
    out.update(Stress=_sym6(S), Cauchy_stress=_sym6(sig),
               VonMises_stress=np.sqrt(1.5 * np.einsum("nij,nij->n", dev, dev)),
               Cauchy_stress_cylindrical=_cylindrical(sig, X + u))
    return out


def _cylindrical(T, pts):
    """Components (rr, tt, zz, rt, tz, zr) of tensors T in the local
    (e_r, e_theta, e_z) basis at the points pts (Cartesian on the axis)."""
    r = np.hypot(pts[:, 0], pts[:, 1])
    axis = r == 0.0
    pts = pts.copy()
    pts[axis, 0], r = 1.0, np.where(axis, 1.0, r)
    Q = np.zeros((len(pts), 3, 3))
    Q[:, 0, :2] = pts[:, :2] / r[:, None]
    Q[:, 1, 0], Q[:, 1, 1] = -pts[:, 1] / r, pts[:, 0] / r
    Q[:, 2, 2] = 1.0
    return _sym6(np.einsum("nik,nkl,njl->nij", Q, T, Q))


def kinematic_fields(X, cells, u, order, fiber):
    """The strain-type struct fields at the nodes (no material model):
    deformation gradient, Jacobian, Green-Lagrange strain (Cartesian and
    cylindrical in the reference basis), fiber stretch and deformed fiber
    direction."""
    F = nodal_deformation_gradient(X, cells, u, order)
    E = 0.5 * (np.einsum("nki,nkj->nij", F, F) - np.eye(3))
    Ff = np.einsum("nij,nj->ni", F, fiber)
    lam = np.linalg.norm(Ff, axis=1)
    return dict(Def_grad=F.reshape(-1, 9), Jacobian=np.linalg.det(F), Strain=_sym6(E),
                Fiber_stretch=lam, Fiber_direction_deformed=Ff / lam[:, None],
                Strain_cylindrical=_cylindrical(E, X))


# through-wall profiles of the struct fields: name -> (field, component)
FIELD_PROFILES = {"J": ("Jacobian", None), "E_RR": ("Strain_cylindrical", 0),
                  "E_TT": ("Strain_cylindrical", 1), "E_ZZ": ("Strain_cylindrical", 2),
                  "sigma_rr": ("Cauchy_stress_cylindrical", 0),
                  "sigma_tt": ("Cauchy_stress_cylindrical", 1),
                  "sigma_zz": ("Cauchy_stress_cylindrical", 2),
                  "sigma_vm": ("VonMises_stress", None), "lambda_f": ("Fiber_stretch", None)}
# the kinematic subset (kinematic_fields)
KINEMATIC_PROFILES = ("J", "E_RR", "E_TT", "E_ZZ", "lambda_f")


def level_means(levels, nodal, keys, use=None):
    """Mean over the nodes of each radial level (only the nodes `use`, a
    boolean mask, if given) of the profile fields keys."""
    use = np.ones(len(levels), bool) if use is None else use
    count = np.bincount(levels[use])
    out = {}
    for key in keys:
        f, c = FIELD_PROFILES[key]
        val = nodal[f] if c is None else nodal[f][:, c]
        out[key] = np.bincount(levels[use], val[use]) / count
    return out


def post(case_dir, fields=True, mid_section=False):
    """Cavity volume and kinematics of a finished case, per load step.

    V is the volume enclosed by the deformed (discrete) endo face. The
    through-wall profiles rho, phi and eta are given at the radial node levels
    R (see `kinematics`), beta is the twist per unit length (zero in a
    symmetric sector). Written to CASE_DIR/kinematics.json. With fields, the
    struct output fields (`struct_fields`) are written per step to
    CASE_DIR/post/fields_<step>.vtu and their through-wall profiles
    (FIELD_PROFILES) added to kinematics.json.

    mid_section compares only the middle element layer (|Z - L/2| <= half an
    element length; the solution of the free-ended tube is symmetric about
    L/2, so its layer means are the mid-length values to O(h^2)): profiles,
    fields and the non-homogeneity from its nodes, beta the local twist rate
    at L/2, and V the mid-section volume pi (Ri + rho(Ri))^2 L of the 0D.
    """
    with open(os.path.join(case_dir, "case.json")) as f:
        case = json.load(f)
    Lm = case.get("mesh_length", case["length"])
    closed = case.get("plate", 0.0) > 0.0
    if closed:
        mesh = mesh_closed_cylinder(case["Ri"], case["Re"], Lm, case["nr"], case["plate"],
                                    case["order"], case["ntheta"], case["nz"], case["nz_plate"],
                                    case["n_in"])
    else:
        mesh = mesh_cylinder(case["Ri"], case["Re"], Lm, case["nr"], case["sector"],
                             case["order"], case["ntheta"], case["nz"])
    X = mesh["points"]
    vol = lambda x: lumen_volume(x, conn, case["sector"]) * case["length"] / Lm  # noqa: E731
    conn = mesh["faces"]["endo"][0]
    dR = (case["Re"] - case["Ri"]) / (case["nr"] * case["order"])
    levels = np.rint((np.hypot(X[:, 0], X[:, 1]) - case["Ri"]) / dR).astype(int)
    files = sorted(glob.glob(os.path.join(case_dir, "results", "result_*.vtu")),
                   key=lambda fn: int(os.path.basename(fn)[7:-4]))  # step number
    if not files:
        raise FileNotFoundError(f"no results in {case_dir}/results")
    section, use = {}, None
    if closed:  # the element layer at the mid-length plane Z = 0
        mid_section = True
        section = dict(z0=0.0, dz=Lm / case["nz"])
    elif mid_section:
        section = dict(z0=0.5 * Lm, dz=0.5 * Lm / case["nz"])
    if mid_section:  # V of the 0D: pi (Ri + rho(Ri))^2 (1 + eps) L
        use = np.abs(X[:, 2] - section["z0"]) <= section["dz"] * (1.0 + 1e-9)
        vol = lambda x, rho_endo, eps=0.0: (  # noqa: E731
            math.pi * (case["Ri"] + rho_endo) ** 2 * (1.0 + eps) * case["length"])
    out = dict(R=(case["Ri"] + dR * np.arange(levels.max() + 1)).tolist(), t=[0.0], P=[0.0],
               V=[vol(X, 0.0) if mid_section else vol(X)], beta=[0.0], eps=[0.0],
               rho_spread=[0.0], section="mid-length" if mid_section else "all")
    zero = [0.0] * (levels.max() + 1)
    out.update(rho=[zero], phi=[zero], eta=[zero])
    if fields:
        from check_materials import svmp_material
        material = svmp_material(os.path.join(case_dir, "solver.xml"))
        m = case["material"]
        fib = cylinder_fibers_at_nodes(X, case["Ri"], case["Re"], m.get("alpha_endo", 0.0),
                                       m.get("alpha_epi", 0.0))
        post_dir = os.path.join(case_dir, "post")
        os.makedirs(post_dir, exist_ok=True)
        u0 = struct_fields(X, mesh["cells"], np.zeros_like(X), case["order"], material,
                           fib["fiber"], fib["sheet"])
        out["fields"] = {k: [v.tolist()]
                         for k, v in level_means(levels, u0, FIELD_PROFILES, use).items()}
    for fn in files:
        n = int(os.path.basename(fn)[7:-4])
        pts, pdata = _read_vtu(fn)
        # results are written in single precision
        if not np.allclose(pts[:, :3], X, rtol=0, atol=1e-6 * case["Re"]):
            raise RuntimeError(f"{fn}: node order differs from the generated mesh")
        x = X + pdata["Displacement"][:, :3]
        t = n * case["dt"]
        k = kinematics(X, x, levels, **section)
        out["t"].append(t)
        out["P"].append(case["P_max"] * min(t / case["T_ramp"], 1.0))
        out["V"].append(vol(x, k["rho"][0], k["eps"] if closed else 0.0) if mid_section
                        else vol(x))
        for key in ("beta", "eps", "rho_spread"):
            out[key].append(k[key])
        for key in ("rho", "phi", "eta"):
            out[key].append(k[key].tolist())
        if fields:
            u = x - X
            sf = struct_fields(X, mesh["cells"], u, case["order"], material, fib["fiber"],
                               fib["sheet"])
            for key, v in level_means(levels, sf, FIELD_PROFILES, use).items():
                out["fields"][key].append(v.tolist())
            pdata = {"Displacement": (u, "Float64", 3)}
            pdata.update({k: (v, "Float64", 1 if v.ndim == 1 else v.shape[1])
                          for k, v in sf.items()})
            write_vtu(os.path.join(post_dir, f"fields_{n:03d}.vtu"), X, mesh["cells"],
                      24 if case["order"] == 2 else 10, pdata, {})
    with open(os.path.join(case_dir, "kinematics.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate", help="write the 3D case of a 0D input")
    g.add_argument("model_0d")
    g.add_argument("case_dir")
    g.add_argument("--nr", type=int, required=True, help="elements through the wall")
    g.add_argument("--formulation", default="struct-p2", choices=sorted(FORMULATIONS))
    g.add_argument("--sector", type=float, choices=[90.0, 360.0],
                   help="default: 90 (symmetry planes), 360 with helical fibers")
    g.add_argument("--half-length", action="store_true",
                   help="mesh half the length (mid-length symmetry plane)")
    g.add_argument("--ntheta", type=int, help="circumferential divisions (default from h)")
    g.add_argument("--nz", type=int, help="axial divisions (default from h)")
    g.add_argument("--steps", type=int, default=20, help="load steps over the ramp")
    g.add_argument("--density", type=float, default=1.0)
    g.add_argument("--ls-type", default="BICG", help="linear solver (BICG, GMRES, ...)")
    g.add_argument("--algebra", default="fsils", help="fsils, trilinos or petsc")
    g.add_argument("--preconditioner", default="fsils", help="e.g. fsils, trilinos-ilut")
    r = sub.add_parser("run", help="run svMultiPhysics on a case")
    r.add_argument("case_dir")
    r.add_argument("--np", type=int, default=1)
    r.add_argument("--exe", help="local svmultiphysics executable instead of Docker")
    r.add_argument("--image", default=IMAGE)
    p = sub.add_parser("post", help="cavity volume and kinematics of a finished case")
    p.add_argument("case_dir")
    a = ap.parse_args()
    if a.cmd == "generate":
        os.makedirs(a.case_dir, exist_ok=True)
        ls = dict(type=a.ls_type, algebra=a.algebra, preconditioner=a.preconditioner)
        c = generate(a.model_0d, a.case_dir, a.nr, a.formulation, a.sector, a.ntheta, a.nz,
                     a.steps, a.density, linear_solver=ls, half_length=a.half_length)
        print(f"{a.case_dir}: {c['n_cells']} tets (order {c['order']}), {c['n_nodes']} nodes, "
              f"nr x ntheta x nz = {c['nr']} x {c['ntheta']} x {c['nz']}")
    elif a.cmd == "run":
        run(a.case_dir, a.np, a.exe, a.image)
    else:
        out = post(a.case_dir)
        for t, P, V, b, rho in zip(out["t"], out["P"], out["V"], out["beta"], out["rho"]):
            print(f"t={t:9.4g}  P={P / MMHG:7.2f} mmHg  V={V * 1e6:10.5f} mL  "
                  f"rho_endo={rho[0] * 1e3:8.5f} mm  rho_epi={rho[-1] * 1e3:8.5f} mm  "
                  f"beta={b:9.5f} rad/m")


if __name__ == "__main__":
    main()
