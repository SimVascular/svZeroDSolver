import os
import json
import numpy as np
import os
import json
import pandas as pd

import sys
sys.path.append(os.path.dirname(__file__))

from .utils import run_test_case_by_name, get_result, run_with_reference, RTOL_FLOW, RTOL_PRES


def test_all():
    '''
    run all test cases and compare against stored reference solution
    '''

    this_file_dir = os.path.abspath(os.path.dirname(__file__))

    results_dir = os.path.join(this_file_dir, 'cases', 'results')

    testfiles = [f for f in os.listdir(os.path.join(this_file_dir, 'cases')) if os.path.isfile(os.path.join(this_file_dir, 'cases', f))]

    # we only want to test the solver, not the calibrator
    testfiles.remove("steadyFlow_calibration.json")

    for file in testfiles:

        with open(os.path.join(this_file_dir, 'cases', file), "r") as f:
            config = json.load(f)

        ref = pd.read_json(os.path.join(results_dir, 'result_' + file))

        run_with_reference(ref, config)
