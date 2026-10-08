# This file contains the stubs for the public Python interface of svZeroDSolver.
# This is necessary for supporting auto-completion in Visual Studio Code or
# PyCharm.
"""svZeroDSolver Python interface."""
from __future__ import annotations
import numpy
import typing
import pandas

__all__ = [
    "Solver",
    "__version__",
    "build_identity",
    "capabilities",
    "calibrate",
    "simulate",
]

__version__: str


def build_identity() -> dict[str, str]:
    """Return the solver version, source revision, and stable build ID."""
    ...


def capabilities() -> dict[str, bool]:
    """Return machine-readable solver capabilities."""
    ...

class Solver:
    """Lumped-parameter solver."""

    @typing.overload
    def __init__(self, arg0: dict) -> None:
        """Create a new lumped-parameter solver.

        Args:
            arg0: Solver configuration dictionary.
        """
        ...
    @typing.overload
    def __init__(self, arg0: str) -> None:
        """Create a new lumped-parameter solver.

        Args:
            arg0: Path to solver configuration file.
        """
        ...
    def get_full_result(self) -> pandas.DataFrame:
        """Get the full result of the simulation.

        Returns:
            Simulation result as a dataframe.
        """
        ...
    def get_single_result(self, arg0: str) -> numpy.ndarray:
        """Get the simulation result for a single degree-of-freedom (DOF).

        Args:
            arg0: Name of the DOF.

        Returns:
            Time-dependent simulation result for DOF.
        """
        ...
    def get_single_result_avg(self, arg0: str) -> float:
        """Get the mean simulation result for a single degree-of-freedom (DOF).

        Args:
            arg0: Name of the DOF.

        Returns:
            Mean simulation result for the DOF.
        """
        ...
    def run(self) -> None:
        """Run the simulation."""
        ...

def calibrate(arg0: dict) -> dict:
    """Run a Levenberg-Marquardt calibration.

    Args:
        arg0: Calibration configuration dictionary.

    Returns:
        Calibrated 0D solver input file.
    """
    ...

@typing.overload
def simulate(arg0: dict) -> pandas.DataFrame:
    """Run a lumped-parameter simulation.

    Args:
        arg0: Simulation configuration file.

    Returns:
            Simulation result as a dataframe.
    """
    ...

@typing.overload
def simulate(arg0: str) -> pandas.DataFrame:
    """Run a lumped-parameter simulation.

    Args:
        arg0: Path to solver configuration file.

    Returns:
            Simulation result as a dataframe.
    """
    ...
