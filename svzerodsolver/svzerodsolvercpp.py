"""Module to call the C++ version of svZeroDSolver.

Expects to find a build in the Release folder.
"""
import importlib.util
import sys
import os
import glob
import orjson
from pandas import read_csv
from io import StringIO

this_file_dir = os.path.abspath(os.path.dirname(__file__))

libfiles = glob.glob(os.path.join(this_file_dir, "..", "Release", "svzerodsolver*.so"))

if not libfiles:
    raise ImportError("No release build of svzerodsolver found.")

spec = importlib.util.spec_from_file_location(
    "svzerodsolvercpp",
    libfiles[0],
)
svzerodsolvercpp = importlib.util.module_from_spec(spec)
sys.modules["svzerodsolvercpp"] = svzerodsolvercpp
spec.loader.exec_module(svzerodsolvercpp)


def run(config):
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
