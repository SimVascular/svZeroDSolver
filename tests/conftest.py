from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import pytest
# from .utils import execute_pysvzerod
import pysvzerod
import json
from os import listdir
from os.path import isfile, join
import os


def pytest_addoption(parser):
    """Pass optional argument to test coverage (slower)"""
    parser.addoption(
        "--coverage",
        action="store_true",
        default="run tests to generate coverage report",
    )


def pytest_configure(config):
    """Workaround because store_true doesn't work"""
    pytest.coverage = config.option.coverage == True


'''
def run_with_reference(
        ref,
        test_config
        ):

    RTOL_PRES = 1.0e-7
    RTOL_FLOW = 1.0e-8

    res = pysvzerod.simulate(test_config)

    if res.shape[1] == 6:
        # we have a result with fields [name, time, p_in, p_out, q_in, q_out]
        for field in ["pressure_in", "pressure_out", "flow_in", "flow_out"]:
            if "pressure" in field:
                assert np.isclose(res[field].to_numpy().all(), ref[field].to_numpy().all(), rtol=RTOL_PRES)
            elif "flow" in field:
                assert np.isclose(res[field].to_numpy().all(), ref[field].to_numpy().all(), rtol=RTOL_FLOW)
    else:
        # we have a result with fields [name, time, y] and the result must be compared based on the name field. name is of format [flow:vessel:outlet]
        for index, row in res.iterrows():
            if "flow" in row["name"]:
                assert np.isclose(row.y, ref.loc[row.name].y, rtol=RTOL_FLOW)
            elif "pressure" in row["name"]:
                assert np.isclose(row.y, ref.loc[row.name].y, rtol=RTOL_PRES)


def test_all():

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


if __name__ == "__main__":
    test_all()

'''