from tempfile import TemporaryDirectory

import numpy as np
import pytest


@pytest.fixture
def tempdir():
    """Temporary directory for test purposes."""
    with TemporaryDirectory() as tempdir:
        yield tempdir


# pass optional argument to test coverage (slower)
def pytest_addoption(parser):
    parser.addoption(
        "--coverage",
        action="store_true",
        default="run tests to generate coverage report",
    )


# weird workaround because store_true doesn't work
def pytest_configure(config):
    pytest.coverage = config.option.coverage == True
