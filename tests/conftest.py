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
