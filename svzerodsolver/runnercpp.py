"""Module to call the C++ version of svZeroDSolver.

Expects to find a build in the Release folder.
"""
import glob
import importlib.util
import os
import sys
from io import StringIO

import orjson
from pandas import read_csv

this_file_dir = os.path.abspath(os.path.dirname(__file__))

libfiles = glob.glob(
    os.path.join(this_file_dir, "..", "Release", "svzerodsolver*.so")
)

if not libfiles:
    raise ImportError(
        "No release build of svzerodsolver found. Please create a build in a folder called 'Release'."
    )

spec = importlib.util.spec_from_file_location(
    "svzerodsolvercpp",
    libfiles[0],
)
svzerodsolvercpp = importlib.util.module_from_spec(spec)
sys.modules["svzerodsolvercpp"] = svzerodsolvercpp
spec.loader.exec_module(svzerodsolvercpp)


def run_from_config(config):
    """Run the C++ svZeroDSolver.

    Args:
        config: Python dict of the configuration.

    Returns:
        Pandas dataframe with the results."""
    result = svzerodsolvercpp.run(
        orjson.dumps(
            config,
            option=orjson.OPT_NAIVE_UTC | orjson.OPT_SERIALIZE_NUMPY,
        )
    )
    return read_csv(StringIO(result))
