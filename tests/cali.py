#!/usr/bin/env python
# coding=utf-8

import pdb
import json
import os
import vtk
import copy
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import svzerodplus
from scipy.interpolate import CubicSpline
from collections import defaultdict
from argparse import ArgumentParser
from io import StringIO

from vtk.util.numpy_support import numpy_to_vtk as n2v
from vtk.util.numpy_support import vtk_to_numpy as v2n

from test_integration_cpp import run_test_case_by_name

FDIR = os.path.dirname(__file__)

# Number of observations to extract in result refinement
NUM_OBS = 100


def refine(x: np.ndarray, y: np.ndarray, num: int) -> np.ndarray:
    """Refine a curve using cubic spline interpolation with derivative.

    Args:
        x: X-coordinates
        y: Y-coordinates
        num: New number of points of the refined data.


    Returns:
        new_y: New y-coordinates
        new_dy: New dy-coordinates
    """
    y = y.copy()
    y[-1] = y[0]
    x_new = np.linspace(x[0], x[-1], num)
    spline = CubicSpline(x, y, bc_type="periodic")
    new_y = spline(x_new)
    new_dy = spline.derivative()(x_new)
    return new_y.tolist(), new_dy.tolist()


def read_geo(fname):
    """
    Read geometry from file, chose corresponding vtk reader
    Args:
        fname: vtp surface or vtu volume mesh

    Returns:
        vtk reader, point data, cell data
    """
    _, ext = os.path.splitext(fname)
    if ext == ".vtp":
        reader = vtk.vtkXMLPolyDataReader()
    elif ext == ".vtu":
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:
        raise ValueError("File extension " + ext + " unknown.")
    reader.SetFileName(fname)
    reader.Update()

    return reader


def collect_arrays(output):
    res = {}
    for i in range(output.GetNumberOfArrays()):
        name = output.GetArrayName(i)
        data = output.GetArray(i)
        res[name] = v2n(data)
    return res


def get_res_names(inp, res_fields):
    # result name list
    res = []

    # get integral for each result
    for i in range(inp.GetPointData().GetNumberOfArrays()):
        res_name = inp.GetPointData().GetArrayName(i)
        field = res_name.split("_")[0]
        num = res_name.split("_")[-1]

        # check if field should be added to output
        if field in res_fields:
            try:
                float(num)
                res += [res_name]
            except ValueError:
                pass

    return res


def load_results_3d(f_res_3d):
    """
    Read 3d results embedded in centerline and sort according to branch at time step
    """
    # read 1d geometry
    reader = read_geo(f_res_3d).GetOutput()
    res = collect_arrays(reader.GetPointData())

    # names of output arrays
    res_names = get_res_names(reader, ["pressure", "velocity"])

    # get time steps
    times = np.unique([float(k.split("_")[1]) for k in res_names])

    # get branch ids
    branches = np.unique(res["BranchId"]).tolist()
    if -1 in branches:
        branches.remove(-1)

    # add time
    out = {"time": times}

    # initilize output arrays [time step, branch]
    for f in res_names:
        name = f.split("_")[0]
        out[name] = {}
        out["path"] = {}
        for br in branches:
            ids = res["BranchId"] == br
            out[name][br] = np.zeros((times.shape[0], np.sum(ids)))
            out["path"][br] = res["Path"][ids]

    # read branch-wise results from geometry
    for f in res_names:
        name, time = f.split("_")
        for br in branches:
            ids = res["BranchId"] == br
            out[name][br][float(time) == times] = res[f][ids]

    # rename velocity to flow
    out["flow"] = out["velocity"]
    del out["velocity"]
    return out


def collect(out, res, name_vessel, name_other, is_out):
    for f in ["flow", "pressure"]:
        if is_out:
            name = ":".join([f, name_vessel, name_other])
            inout = "_out"
        else:
            name = ":".join([f, name_other, name_vessel])
            inout = "_in"
        i_vesel = res["name"] == name_vessel
        data = np.array(res[f + inout][i_vesel])
        time = np.array(res["time"][i_vesel])
        out["y"][name], out["dy"][name] = refine(time, data, NUM_OBS)


def collect_results_from_0d(inp, res):
    vessels = np.array(inp["vessels"])
    out = {"y": {}, "dy": {}}
    # collect connections between branches and junctions
    for jc in inp["junctions"]:
        for io in ["inlet", "outlet"]:
            for vs in vessels[jc[io + "_vessels"]]:
                collect(out, res, vs["vessel_name"], jc["junction_name"], io == "inlet")
    # collect connections between branches and boundary conditions
    for vs in inp["vessels"]:
        if "boundary_conditions" in vs:
            for bc_type, bc_name in vs["boundary_conditions"].items():
                collect(out, res, vs["vessel_name"], bc_name, bc_type == "outlet")
    
    # export results, time, and flow on same resampled time discretization
    flow = out["y"]["flow:INFLOW:branch0_seg0"]
    return out, np.linspace(0, np.max(res["time"]), len(flow)).tolist(), flow


def collect_results_from_3d(inp, res):
    # scale time to 0d
    for bc in inp["boundary_conditions"]:
        if bc["bc_name"] == "INFLOW":
            t_max = bc["bc_values"]["t"][-1]
    # todo: pick only last cycle (still doesn't work when doing it manually)
    nt = 20
    time = np.linspace(0, t_max, nt)

    # convert 3d results to 0d result format
    res0d = defaultdict(list)
    for vessel in inp["vessels"]:
        res0d["name"] += [vessel["vessel_name"]] * nt
        res0d["time"] += time.tolist()

        # get branch and segment id
        names = vessel["vessel_name"].split("_")
        br, seg = [int(names[i].split(n)[-1]) for i, n in enumerate(["branch", "seg"])]

        # extract results
        if seg == 0:
            # reset branch length
            br_len = 0.0
        for io in ["_in", "_out"]:
            # find closest centerline point
            loc = np.argmin(np.abs(res["path"][br] - br_len))
            for f in ["flow", "pressure"]:
                res0d[f + io] += res[f][br][-nt:, loc].tolist()
                # pdb.set_trace()
            br_len += vessel["vessel_length"]

    # convert to numpy arrays
    for k, v in res0d.items():
        res0d[k] = np.array(v)
    return collect_results_from_0d(inp, res0d)


def compare(geo, dim):
    print("\n", geo, "\n\n")

    # run the estimation
    inp0, cali = estimate(geo, dim)

    # compare 0D element values
    for i in range(len(inp0["vessels"])):
        for ele in cali["vessels"][i]["zero_d_element_values"].keys():
            sol = cali["vessels"][i]["zero_d_element_values"][ele]
            if ele in inp0["vessels"][i]["zero_d_element_values"]:
                ref = inp0["vessels"][i]["zero_d_element_values"][ele]
            else:
                ref = 0.0

            # calculate error
            if ref == 0.0:
                err = 1.0
            else:
                err = np.abs(sol / ref - 1.0)

            # print results
            out = str(i) + "\t" + ele[0]
            for j in [ref, sol]:
                out += "\t\t{:.1e}".format(j)
            out += "\t\t{:+d}".format(int(np.log(err)))
            print(out)


def plot(dim):
    files = np.loadtxt("geometries_paper.txt", dtype="str")
    assert len(files) == 72, "wrong number of files"

    # compare 0D element values
    nx = 8
    ny = 9
    colors = {"C": "b", "L": "g", "R": "r", "s": "k"}
    fig, ax = plt.subplots(nx, ny, figsize=(ny * 2, nx * 2), dpi=300)

    for j, fname in enumerate(files):
        print(fname)
        ab = np.unravel_index(j, (nx, ny))
        # read results
        end = "_from_" + str(dim) + "d.json"
        with open(os.path.join(FDIR, "cases", "vmr", fname + ".json")) as f:
            inp = json.load(f)
        with open(os.path.join(FDIR, "cases", "vmr", fname + "_optimal" + end)) as f:
            opt = json.load(f)

        for i in range(len(inp["vessels"])):
            for ele in opt["vessels"][i]["zero_d_element_values"].keys():
                sol = opt["vessels"][i]["zero_d_element_values"][ele]
                if ele in inp["vessels"][i]["zero_d_element_values"]:
                    ref = inp["vessels"][i]["zero_d_element_values"][ele]
                else:
                    ref = 0.0
                ax[ab].loglog(ref, sol, "o", color=colors[ele[0]])
        ax[ab].set_title(fname)
        ax[ab].set_aspect("equal", adjustable="datalim")
        ax[ab].set_xticklabels([])
        ax[ab].set_yticklabels([])
        ax[ab].grid(True)
        xlim = ax[ab].get_xlim()
        ylim = ax[ab].get_ylim()
        lim = []
        for i in range(2):
            if xlim[i] > ylim[i]:
                lim += [xlim[i]]
            else:
                lim += [ylim[i]]
        ax[ab].plot(lim, lim, color="k")
        ax[ab].set_xlim(xlim)
        ax[ab].set_ylim(ylim)
    plt.tight_layout()
    fig.savefig("calibration_" + str(dim) + "d.png", bbox_inches="tight")


def estimate(geo, dim):
    # read input file
    with open(os.path.join(FDIR, "cases", "vmr", geo + ".json")) as f:
        inp = json.load(f)
    inp0 = copy.deepcopy(inp)

    # get forward results
    if dim == 0:
        # run simulation
        res = run_test_case_by_name(geo, output_variable_based=True, folder="vmr")

        # collect results in format required for calibrator
        out, time, flow = collect_results_from_0d(inp, res)
    elif dim == 3:
        # read 3d results from centerline
        res = load_results_3d(os.path.join("cases", "vmr_from_3d", geo + ".vtp"))
        out, time, flow = collect_results_from_3d(inp, res)
    else:
        raise ValueError("unknown dimension: " + str(dim))
    
    # replace inflow to match calibraiton data
    # todo: this is probably not necessary
    for bc in inp["boundary_conditions"]:
        if bc["bc_name"] == "INFLOW":
            bc["bc_values"]["t"] = time
            bc["bc_values"]["Q"] = flow

    # set all elements to zero
    for vessel in inp["vessels"]:
        for ele in vessel["zero_d_element_values"].keys():
            vessel["zero_d_element_values"][ele] = 0.0

    # add calibration parameters
    inp["calibration_parameters"] = {
        "tolerance_gradient": 1e-5,
        "tolerance_increment": 1e-9,
        "maximum_iterations": 20,
        "calibrate_stenosis_coefficient": True,
        "set_capacitance_to_zero": False,
    }

    # only calibrate to last cycle
    inp["simulation_parameters"]["output_all_cycles"] = False

    # write output json
    end = "_from_" + str(dim) + "d"
    fname_out = os.path.join(FDIR, "cases", "vmr", geo + "_calibrate" + end + ".json")
    with open(fname_out, "w") as f:
        json.dump(inp | out, f, indent=4)

    # run the calibration
    with open(fname_out) as ff:
        config = json.load(ff)
    try:
        cali = svzerodplus.calibrate(config)
    except RuntimeError as e:
        print("calibration failed: ", e)
        return

    # write calibrated json
    fname_out = os.path.join(FDIR, "cases", "vmr", geo + "_optimal" + end)
    with open(fname_out + ".json", "w") as f:
        json.dump(cali, f, indent=4)

    # run 0d simulation with optimal parameters
    # todo: compare res and res2
    # res2 = run_test_case_by_name(fname_out, output_variable_based=False, folder="vmr")
    return inp0, cali


if __name__ == "__main__":
    p = ArgumentParser(description="Run inverse problem from forward test case")
    p.add_argument("dim", type=int, choices=[0, 3], help="input simulation dimension")
    p.add_argument("case", nargs="?", default=None, help="test case (loop all if none)")
    args = p.parse_args()

    if args.case:
        # run user-provided test case
        compare(args.case, args.dim)
    else:
        # loop over all test cases
        files = os.path.join(FDIR, "cases", "vmr", "*.json")
        sample = sorted(glob.glob(files))
        for case in sample:
            if "calibrate" not in case and "optimal" not in case:
                compare(os.path.basename(os.path.splitext(case)[0]), args.dim)
        plot(args.dim)
