import json
import os
import subprocess
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

# global boolean to perform coverage testing
# (run executables instead of Python interface, much slower)
from pytest import coverage
# coverage = False

import pysvzerod

this_file_dir = os.path.abspath(os.path.dirname(__file__))

RTOL_PRES = 1.0e-7
RTOL_FLOW = 1.0e-7

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

def compare_result_with_reference(res, ref, rtol_pres=1.0e-7, rtol_flow=1.0e-7, output_variable_based=False):
    '''
    Compare the result with the reference.

    Args:
        res: result as pandas DataFrame 
        ref: reference result as pandas DataFrame
        output_variable_based: whether to compare based on columns (True) or row-based with 'name' field (False)

    Returns:
        pd.DataFrame with columns: ["variable", "name", "y_expected", "y_actual", "rel_diff - tol", "within_tolerance"]
    '''

    results = []

    if output_variable_based:
        assert len(res) == len(ref), "Result and reference must have the same number of rows"

        name = ref["name"]
        y_expected = ref["y"]
        t_reference = ref["time"]
        y_actual = res["y"]
        t_result = res["time"]
        tol = name.map(lambda n: rtol_flow if "flow" in n else rtol_pres)

        diff_vs_zero = abs(y_actual - y_expected) - tol - (tol * abs(y_expected))
        within_tol = diff_vs_zero <= 0.0

        for var_name, t_e, y_e, t_a, y_a, rd, wt in zip(name, t_reference, y_expected, t_result, y_actual, diff_vs_zero, within_tol):
            variable, block_name = var_name.split(":", 1)
            results.append({
                "variable": variable,
                "name": block_name,
                "t_reference": t_e,
                "y_expected": y_e,
                "t_result": t_a,
                "y_result": y_a,
                "rel_diff - tol": rd,
                "within_tolerance": wt
            })

    else:
        for col in res.columns:
            tol = rtol_pres if "pressure" in col else rtol_flow if "flow" in col else None
            if tol is None:
                    continue

            y_expected = ref[col]
            y_actual = res[col]
            diff_vs_zero = abs(y_actual - y_expected) - tol - (tol * abs(y_expected))
            within_tol = diff_vs_zero <= 0.0

            for idx, y_e, y_a, rd, wt in zip(res.index, y_expected, y_actual, diff_vs_zero, within_tol):
                name = res.loc[idx, "name"] if "name" in res.columns else str(idx)
                t_e = ref.loc[idx, "time"] if "time" in ref.columns else None
                t_a = res.loc[idx, "time"] if "time" in res.columns else None
                results.append({
                    "variable": col,
                    "name": name,
                    "t_reference": t_e,
                    "y_expected": y_e,
                    "t_result": t_a,
                    "y_result": y_a,
                    "rel_diff - tol": rd,
                    "within_tolerance": wt
                })

    return pd.DataFrame(results, columns=[
        "variable", "name", "t_reference", "y_expected", "t_result", "y_result", "rel_diff - tol", "within_tolerance"
    ])


def run_with_reference(
        ref,
        test_config,
        rtol_pres=1.0e-7,
        rtol_flow=1.0e-7,
        ):


    res, config = execute_pysvzerod(test_config, "solver")

    output_variable_based = config["simulation_parameters"].get("output_variable_based", False)

    difference = compare_result_with_reference(res, ref, rtol_pres, rtol_flow, output_variable_based)

    if not difference["within_tolerance"].all():
        # Extract only differing rows for a cleaner error message
        differing_rows = difference[~difference["within_tolerance"]]
        if not differing_rows.empty:
            print("Test failed in the following rows:\n", differing_rows.to_string(index=False))
            raise AssertionError("Differences exceed tolerance.")
        else:
            raise AssertionError("Differences exceed tolerance but no specific rows found.")



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


# from scipy.interpolate import interp1d

# def resample_external_solver_inputs():
#     """
#     Load JSON files, resample external_solver_coupling_blocks["values"]["t"] and ["Q"]
#     to a uniform time grid defined by simulation_parameters["number_of_time_pts_per_cardiac_cycle"]
#     and number_of_cardiac_cycles.

#     """
#     cardiac_cycle_period = 0.99

#     json_paths = [os.path.join(this_file_dir, "cases", "coupledBlock_closedLoopHeart_singleVessel.json"),
#                   os.path.join(this_file_dir, "cases", "coupledBlock_closedLoopHeart_withCoronaries.json"),]

#     updated_jsons = {}

#     for path in json_paths:
#         with open(path, "r") as f:
#             data = json.load(f)

#         sim_params = data.get("simulation_parameters", {})
#         num_pts_per_cycle = sim_params.get("number_of_time_pts_per_cardiac_cycle", 100)
#         total_pts = num_pts_per_cycle

#         # set cardiac cycle period
#         data["closed_loop_blocks"][0]["cardiac_cycle_period"] = cardiac_cycle_period

#         # Estimate total time from last point in original t
#         for block in data.get("external_solver_coupling_blocks", []):
#             t_original = np.array(block["values"]["t"], dtype=np.float64)
#             q_original = np.array(block["values"]["Q"], dtype=np.float64)

#             # Create a uniform time grid from t[0] to t[-1]
#             t_uniform = np.linspace(0.0, cardiac_cycle_period, total_pts, endpoint=True)
#             t_uniform = np.floor(t_uniform * 1000) / 1000

#             # Interpolate Q to match new time points
#             interp_q = interp1d(t_uniform, q_original, kind='linear', fill_value="extrapolate")
#             q_uniform = interp_q(t_uniform)

#             # Update the JSON block
#             block["values"]["t"] = t_uniform.tolist()
#             block["values"]["Q"] = q_uniform.tolist()

#         updated_jsons[path] = data
    
#     # Save the updated JSON files
#     for path, updated_data in updated_jsons.items():
#         with open(path, "w") as f:
#             json.dump(updated_data, f, indent=4)
#             print(f"Updated JSON saved to {path}")

# if __name__ == "__main__":

#     resample_external_solver_inputs()