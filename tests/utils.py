import json
import os
import subprocess
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

# global boolean to perform coverage testing
# (run executables instead of Python interface, much slower)
from pytest import coverage

import svzerodplus

this_file_dir = os.path.abspath(os.path.dirname(__file__))

RTOL_PRES = 1.0e-7
RTOL_FLOW = 1.0e-8


def execute_svzerodplus(testfile, mode):
    """Execute svzerodplus (via Python interface or executable).

    Args:
        testfile: Path to the input file.
        mode: svZeroDPlus application (solver or calibrator).
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
            result = svzerodplus.simulate(config)
        elif mode == "calibrator":
            result = svzerodplus.calibrate(config)

    return result, config


def run_test_case_by_name(name, output_variable_based=False, folder="."):
    """Run a test case by its case name.

    Args:
        name: Name of the test case.
        testdir: Directory for performing the simulation.
    """
    # file name of test case
    testfile = os.path.join(this_file_dir, "cases", name + ".json")

    # run test
    result, config = execute_svzerodplus(testfile, "solver")

    if output_variable_based == False:
        output = {
            "pressure_in": {},
            "pressure_out": {},
            "flow_in": {},
            "flow_out": {},
        }

        last_seg_id = 0

        for vessel in config["vessels"]:
            name = vessel["vessel_name"]
            branch_id, seg_id = name.split("_")
            branch_id, seg_id = int(branch_id[6:]), int(seg_id[3:])
            vessel_id = vessel["vessel_id"]

            if seg_id == 0:
                output["pressure_in"][branch_id] = np.array(
                    result[result.name == name]["pressure_in"]
                )
                output["flow_in"][branch_id] = np.array(
                    result[result.name == name]["flow_in"]
                )
                output["pressure_out"][branch_id] = np.array(
                    result[result.name == name]["pressure_out"]
                )
                output["flow_out"][branch_id] = np.array(
                    result[result.name == name]["flow_out"]
                )
            elif seg_id > last_seg_id:
                output["pressure_out"][branch_id] = np.array(
                    result[result.name == name]["pressure_out"]
                )
                output["flow_out"][branch_id] = np.array(
                    result[result.name == name]["flow_out"]
                )

            last_seg_id = seg_id

    elif output_variable_based == True:
        output = result

    return output


def get_result(result_array, field, branch, time_step):
    """ "Get results at specific field, branch, branch_node and time step."""
    # extract result
    return result_array[field][branch][time_step]