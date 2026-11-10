# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Myofiber directions of a thick-walled cylinder for svMultiPhysics.

Writes fiber.vtu, sheet.vtu and normal.vtu next to the mesh: copies of the
mesh with the direction in the array FIB_DIR, per node (point data, read with
Fiber_direction_from_file_points) or per element (cell data,
Fiber_direction_from_file_cells), in the basis convention of svMultiPhysics'
rule-based fiber generator (utilities/fiber_generation, Bayer et al. 2012):

  e_t  transmural, endo -> epi            = e_R
  e_l  longitudinal                       = e_Z
  e_c  circumferential, e_c = e_l x e_t   = e_Theta
  f = cos(a) e_c + sin(a) e_l,   n = -sin(a) e_c + cos(a) e_l,   s = e_t

with the helix angle a rotating about e_t (no transverse angle). As in the
ChamberCylinder 0D block, a is linear in the radius,
a(R) = a_endo (Re - R)/(Re - Ri) + a_epi (R - Ri)/(Re - Ri),
evaluated at the nodes or at the element centroids. (The generator's Laplace
coordinate would be ln(R/Ri)/ln(Re/Ri) for a cylinder, a different profile.)
svMultiPhysics interpolates nodal directions to the quadrature points (principal
direction of the interpolated dyad), so nodal data follow the curved basis
vectors to the order of the element; cell data are constant per element.

Usage: python cylinder_fibers.py MESH.vtu --model-0d MODEL_0D.json
       python cylinder_fibers.py MESH.vtu --Ri 0.019 --Re 0.033 --alpha-endo 60 --alpha-epi -60
"""
import argparse
import json
import math
import os

import numpy as np


def cylinder_fibers(points, cells, Ri, Re, alpha_endo, alpha_epi, location="points"):
    """Fiber, sheet and sheet-normal unit vectors, each (N, 3), at the nodes
    (location "points") or at the cell centroids of the four corner nodes
    ("cells", also for 10-node cells)."""
    c = points if location == "points" else points[cells[:, :4]].mean(axis=1)
    R = np.hypot(c[:, 0], c[:, 1])
    t = (R - Ri) / (Re - Ri)
    a = np.radians(alpha_endo * (1.0 - t) + alpha_epi * t)
    eR = np.stack([c[:, 0] / R, c[:, 1] / R, np.zeros_like(R)], axis=1)
    eT = np.stack([-eR[:, 1], eR[:, 0], np.zeros_like(R)], axis=1)
    eZ = np.tile([0.0, 0.0, 1.0], (len(R), 1))
    ca, sa = np.cos(a)[:, None], np.sin(a)[:, None]
    return dict(fiber=ca * eT + sa * eZ, sheet=eR, normal=-sa * eT + ca * eZ)


def write_fibers(folder, points, cells, vtk_type, fibers, location="points"):
    """fiber.vtu, sheet.vtu, normal.vtu with FIB_DIR as point or cell data;
    returns the paths."""
    from cylinder_3d import write_vtu
    paths = []
    for name in ("fiber", "sheet", "normal"):
        path = os.path.join(folder, f"{name}.vtu")
        data = {"FIB_DIR": (fibers[name], "Float64", 3)}
        if location == "points":
            write_vtu(path, points, cells, vtk_type, data, {})
        else:
            write_vtu(path, points, cells, vtk_type, {}, data)
        paths.append(path)
    return paths


def _read_mesh(path):
    from vtkmodules.util.numpy_support import vtk_to_numpy
    from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
    r = vtkXMLUnstructuredGridReader()
    r.SetFileName(path)
    r.Update()
    g = r.GetOutput()
    types = vtk_to_numpy(g.GetCellTypesArray())
    if len(set(types)) != 1:
        raise ValueError(f"{path}: mixed cell types")
    conn = vtk_to_numpy(g.GetCells().GetConnectivityArray())
    return vtk_to_numpy(g.GetPoints().GetData()), conn.reshape(g.GetNumberOfCells(), -1), \
        int(types[0])


def main():
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("mesh", help="volume mesh (.vtu)")
    ap.add_argument("--model-0d", help="0D input with the ChamberCylinder geometry and angles")
    ap.add_argument("--Ri", type=float)
    ap.add_argument("--Re", type=float)
    ap.add_argument("--alpha-endo", type=float)
    ap.add_argument("--alpha-epi", type=float)
    ap.add_argument("--location", choices=["points", "cells"], default="points")
    a = ap.parse_args()
    if a.model_0d:
        with open(a.model_0d) as f:
            cfg = json.load(f)
        v = next(b for b in cfg["vessels"] if b["zero_d_element_type"] == "ChamberCylinder")
        v = v["zero_d_element_values"]
        Ri, Re, ae, ap_ = v["Ri"], v["Re"], v["alpha_endo"], v["alpha_epi"]
    else:
        Ri, Re, ae, ap_ = a.Ri, a.Re, a.alpha_endo, a.alpha_epi
    pts, cells, vtk_type = _read_mesh(a.mesh)
    fib = cylinder_fibers(pts, cells, Ri, Re, ae, ap_, a.location)
    for p in write_fibers(os.path.dirname(os.path.abspath(a.mesh)), pts, cells, vtk_type, fib,
                          a.location):
        print(p)
    a0 = math.degrees(math.asin(np.clip(fib["fiber"][:, 2], -1, 1).max()))
    print(f"{len(cells)} cells, {a.location} data, helix angle up to {a0:.1f} deg")


if __name__ == "__main__":
    main()
