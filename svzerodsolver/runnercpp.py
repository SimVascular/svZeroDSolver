# Copyright (c) Stanford University, The Regents of the University of
#               California, and others.
#
# All Rights Reserved.
#
# See Copyright-SimVascular.txt for additional details.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""Module to call the C++ version of svZeroDSolver.

Expects to find a build in the Release folder.
"""
import glob
import importlib.util
import os
import sys
from io import BytesIO

import orjson
from pandas import read_csv

this_file_dir = os.path.abspath(os.path.dirname(__file__))

libfiles = glob.glob(
    os.path.join(this_file_dir, "..", "Release", "libsvzerodsolver*.so")
)

if not libfiles:
    raise ImportError(
        "No release build of svzerodsolver found. Please create a build in a folder called 'Release'."
    )

spec = importlib.util.spec_from_file_location(
    "libsvzerodsolver",
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
    return read_csv(BytesIO(result))
