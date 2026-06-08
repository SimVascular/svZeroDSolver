"""
Shared simulation runner for sv0D Tuning Framework.

Runs sv0D with given parameter values and returns the solver and output extractor.
Used by both SV0DTuner and SensitivityAnalyzer to avoid duplicating simulation logic.
"""

import numpy as np
import pysvzerod
from typing import List, Dict, Any, Tuple

from .parameter_handler import ParameterHandler
from .output_extractor import OutputExtractor


def run_simulation(
    param_handler: ParameterHandler,
    parameters: List[Dict[str, Any]],
    param_values: np.ndarray,
) -> Tuple[pysvzerod.Solver, OutputExtractor]:
    """
    Run sv0D simulation with given parameter values.

    Args:
        param_handler: ParameterHandler instance with model config.
        parameters: List of parameter config dicts with 'name' (and optionally 'bounds').
        param_values: Array of parameter values in same order as parameters.

    Returns:
        Tuple of (solver, extractor). Solver has run() already called.
        Extractor is built from the solver.

    Raises:
        Exception: If parameter update or simulation fails.
    """
    param_names = [p["name"] for p in parameters]
    for name, value in zip(param_names, param_values):
        if isinstance(value, np.ndarray):
            value = float(value.item())
        elif not isinstance(value, (int, float)):
            value = float(value)
        param_handler.set_parameter(name, value)

    config_dict = param_handler.get_config()
    solver = pysvzerod.Solver(config_dict)
    solver.run()
    extractor = OutputExtractor(solver)
    return solver, extractor
