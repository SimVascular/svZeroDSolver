from tempfile import TemporaryDirectory

import numpy as np
import pytest
from .utils import execute_pysvzerod


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


def run_with_reference(
        ref,
        test_config
        ):

    RTOL_PRES = 1.0e-7
    RTOL_FLOW = 1.0e-8
    
    res, config = execute_pysvzerod(test_config, "solver")

    for field in ["pressure_in", "pressure_out", "flow_in", "flow_out"]:
        if "pressure" in field:
            assert np.isclose(res[field].to_numpy(), ref[field].to_numpy(), rtol=RTOL_PRES)
        elif "flow" in field:
            assert np.isclose(res[field].to_numpy(), ref[field].to_numpy(), rtol=RTOL_FLOW)

    print('test passed!')