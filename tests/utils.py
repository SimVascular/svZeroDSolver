import json
import os
import subprocess
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

# global boolean to perform coverage testing
# (run executables instead of Python interface, much slower)
from pytest import coverage

import pysvzerod

this_file_dir = os.path.abspath(os.path.dirname(__file__))

RTOL_PRES = 1.0e-7
RTOL_FLOW = 1.0e-8


def execute_pysvzerod(testfile, mode):
    """Execute pysvzerod (via Python interface or executable).

    Args:
        testfile: Path to the input file.
        mode: svZeroDSolver application (solver or calibrator).
    """
    assert mode in ["solver", "calibrator"], "unknown mode: " + mode

    # read configuration
    with open(testfile) as ff:
        config = json.load(ff)

    if coverage:
        # run via executable (slow)
        with TemporaryDirectory() as tempdir:
            out_name = os.path.join(tempdir, "out")
            exe = os.path.join(this_file_dir, "..", "Release", "svzerod")
            subprocess.run([exe + mode, testfile, out_name])
            if mode == "solver":
                result = pd.read_csv(out_name)
            elif mode == "calibrator":
                with open(out_name) as ff:
                    result = json.load(ff)
    else:
        # run via Python binding (fast)
        if mode == "solver":
            result = pysvzerod.simulate(config)
        elif mode == "calibrator":
            result = pysvzerod.calibrate(config)

    return result, config


def run_with_reference(
        ref,
        test_config
        ):


    res, config = execute_pysvzerod(test_config, "solver")

    if res.shape[1] == 6:
        # we have a result with fields [name, time, p_in, p_out, q_in, q_out]
        for field in ["pressure_in", "pressure_out", "flow_in", "flow_out"]:
            if "pressure" in field:
                assert np.isclose(res[field].to_numpy().all(), ref[field].to_numpy().all(), rtol=RTOL_PRES)
            elif "flow" in field:
                assert np.isclose(res[field].to_numpy().all(), ref[field].to_numpy().all(), rtol=RTOL_FLOW)
    else:
        # we have a result with fields [name, time, y] and the result must be compared based on the name field. name is of format [flow:vessel:outlet]
        # we will compare the average of each branch
        avg_res_flow = []
        avg_ref_flow = []
        avg_res_pres = []
        avg_ref_pres = []
        for index, row in res.iterrows():

            if "flow" in row["name"]:
                if row["name"] == res.iloc[index + 1]["name"]:
                    # we are compilng the results for a branch
                    avg_ref_flow.append(ref.loc[row.name].y)
                    avg_res_flow.append(row.y)
                elif avg_res_flow == []:
                    # there is only one result for this branch
                    assert np.isclose(row.y, ref.loc[row.name].y, rtol=RTOL_FLOW)
                else:
                    # we are on the last result for this branch
                    avg_ref_flow.append(ref.loc[row.name].y)
                    avg_res_flow.append(row.y)
                    assert np.isclose(np.array(avg_res_flow).all(), np.array(avg_ref_flow).all(), rtol=RTOL_FLOW)
                    avg_res_flow = []
                    avg_ref_flow = []
                    
            elif "pressure" in row["name"]:
                if index == len(res) - 1:
                    # we are on the last row
                    avg_ref_pres.append(ref.loc[row.name].y)
                    avg_res_pres.append(row.y)
                    assert np.isclose(np.array(avg_res_pres).all(), np.array(avg_ref_pres).all(), rtol=RTOL_PRES)
                elif row["name"] == res.iloc[index + 1]["name"]:
                    # we are compilng the results for a branch
                    avg_ref_pres.append(ref.loc[row.name].y)
                    avg_res_pres.append(row.y)
                elif avg_res_pres == []:
                    # there is only one result for this branch
                    assert np.isclose(row.y, ref.loc[row.name].y, rtol=RTOL_PRES)
                else:
                    # we are on the last result for this branch
                    avg_ref_pres.append(ref.loc[row.name].y)
                    avg_res_pres.append(row.y)
                    # round the result to 10 decimal places to avoid floating point errors with reference solution computed on ubuntu OS
                    assert np.isclose(np.array(avg_res_pres).round(10).all(), np.array(avg_ref_pres).all(), rtol=RTOL_PRES)
                    avg_res_pres = []
                    avg_ref_pres = []


def run_test_case_by_name(name, output_variable_based=False, folder="."):
    """Run a test case by its case name.

    Args:
        name: Name of the test case.
        testdir: Directory for performing the simulation.
    """
    # file name of test case
    testfile = os.path.join(this_file_dir, "cases", name + ".json")

    # run test
    result, config = execute_pysvzerod(testfile, "solver")

    if not output_variable_based:
        output = {
            "pressure_in": {},
            "pressure_out": {},
            "flow_in": {},
            "flow_out": {},
        }

        for vessel in config["vessels"]:
            name = vessel["vessel_name"]
            branch_id, seg_id = name.split("_")
            branch_id, seg_id = int(branch_id[6:]), int(seg_id[3:])

            for f in ["pressure", "flow"]:
                for l in ["in"] * (seg_id == 0) + ["out"]:
                    n = f + "_" + l
                    ids = result.name == name
                    output[n][branch_id] = np.array(result[ids][n])

    else:
        output = result

    return output


def get_result(result_array, field, branch, time_step):
    """ "Get results at specific field, branch, branch_node and time step."""
    # extract result
    return result_array[field][branch][time_step]
