#!/usr/bin/env python
# coding=utf-8

import pdb
import json
import os
import copy
import glob
import argparse
import numpy as np
import svzerodplus

from test_integration_cpp import run_test_case_by_name

def collect_results(inp, res):
    out = {"y": {}, "dy": {}}
    for f in ["flow", "pressure"]:
        # loop over all junctions
        for junction in inp["junctions"]:
            # junction inlets
            for br in junction["inlet_vessels"]:
                for vessel in inp["vessels"]:
                    if vessel["vessel_id"] == br:
                        name = ":".join([f, vessel["vessel_name"], junction["junction_name"]])
                        m = list(res[f + "_out"][res["name"] == vessel["vessel_name"]])
                        out["y"][name] = m
                        out["dy"][name] = np.gradient(m).tolist()
            # junction outlets
            for br in junction["outlet_vessels"]:
                for vessel in inp["vessels"]:
                    if vessel["vessel_id"] == br:
                        name = ":".join([f, junction["junction_name"], vessel["vessel_name"]])
                        m = list(res[f + "_in"][res["name"] == vessel["vessel_name"]])
                        out["y"][name] = m
                        out["dy"][name] = np.gradient(m).tolist()
        # loop over all boundary conditions
        for vessel in inp["vessels"]:
            if "boundary_conditions" in vessel:
                if "inlet" in vessel["boundary_conditions"]:
                    name = ":".join([f, vessel["boundary_conditions"]["inlet"], vessel["vessel_name"]])
                    m = list(res[f + "_in"][res["name"] == vessel["vessel_name"]])
                    out["y"][name] = m
                    out["dy"][name] = np.gradient(m).tolist()
                if "outlet" in vessel["boundary_conditions"]:
                    name = ":".join([f, vessel["vessel_name"], vessel["boundary_conditions"]["outlet"]])
                    m = list(res[f + "_out"][res["name"] == vessel["vessel_name"]])
                    out["y"][name] = m
                    out["dy"][name] = np.gradient(m).tolist()
    return out

def compare(fname):
    print("\n", fname, "\n")

    # read input file
    with open(os.path.join("cases", fname + '.json')) as f:
        inp = json.load(f)
    inp0 = copy.deepcopy(inp)

    # set all elements to zero
    for vessel in inp["vessels"]:
        for ele in vessel['zero_d_element_values'].keys():
            vessel['zero_d_element_values'][ele] = 0.0

    # run simulation
    res = run_test_case_by_name(fname, output_variable_based=True)

    # collect results in format required for calibrator
    out = collect_results(inp, res)

    # add calibration parameters
    inp["calibration_parameters"] = {"tolerance_gradient": 1e-6,
                                     "tolerance_increment": 1e-9,
                                     "maximum_iterations": 200,
                                     "calibrate_stenosis_coefficient": True,
                                     "set_capacitance_to_zero": False}

    # write output json
    fname_out = os.path.join("cases", fname + "_calibration.json")
    with open(fname_out, 'w') as f:
        json.dump(inp | out, f, indent=4)
    
    # run the calibration
    with open(fname_out) as ff:
        config = json.load(ff)
    try:
        cali = svzerodplus.calibrate(config)
        os.remove(fname_out)
    except RuntimeError as e:
        print("calibration failed: ", e)
        os.remove(fname_out)
        return

    # compare 0D element values
    for i in range(len(inp0["vessels"])):
        for ele in cali["vessels"][i]["zero_d_element_values"].keys():
            if ele in inp0["vessels"][i]["zero_d_element_values"]:
                ref = inp0["vessels"][i]["zero_d_element_values"][ele]
            else:
                ref = 0.0
            print(i, ele, ref, cali["vessels"][i]["zero_d_element_values"][ele])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run an inverse problem from a forward test case")
    parser.add_argument("case", nargs='?', default=None, help="test case name (loop all if none provided)")
    args = parser.parse_args()

    if args.case:
        # run user-provided test case
        compare(args.case)
    else:
        # loop over all test cases
        for fn in sorted(glob.glob(os.path.join("cases", "*.json"))):
            if "calibration" not in fn:
                compare(os.path.basename(os.path.splitext(fn)[0]))