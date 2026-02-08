"""
Objective function for sv0D Tuning Framework.
"""

import numpy as np
from typing import Dict, List, Callable, Optional, Union
from scipy.interpolate import interp1d


def _interpolate_to_target_times(
    sim_times: np.ndarray,
    sim_values: np.ndarray,
    target_times: np.ndarray
) -> np.ndarray:
    """Interpolate simulated values to target time points."""
    valid_mask = np.isfinite(sim_values) & np.isfinite(sim_times)
    if np.sum(valid_mask) < 2:
        return np.full_like(target_times, np.nan)
    interp_func = interp1d(
        sim_times[valid_mask],
        sim_values[valid_mask],
        kind='linear',
        bounds_error=False,
        fill_value='extrapolate'
    )
    return interp_func(target_times)


def _l2_error(
    target_values: np.ndarray,
    sim_values: np.ndarray,
    normalize: bool
) -> float:
    """
    L2 error between target and simulated values. Works for both time series
    (length N) and scalars (length 1).
    """
    target_values = np.asarray(target_values)
    sim_values = np.asarray(sim_values)
    error = float(np.linalg.norm(sim_values - target_values))
    if normalize:
        norm = np.linalg.norm(target_values)
        if norm > 1e-14:
            return error / norm
        return error
    return error


class ObjectiveFunction:
    """
    Objective function for optimization. Computes weighted error between
    simulated values and targets. Supports scalars and time series.
    """

    def __init__(
        self,
        targets: List[Dict],
        normalize: bool = False,
        custom_function: Optional[Callable] = None
    ):
        """
        Initialize objective function.

        Args:
            targets: List of target specifications
            normalize: If True, use relative error (normalized by target)
            custom_function: Optional custom function (overrides default)
        """
        self.targets = targets
        self.normalize = normalize
        self.custom_function = custom_function
        self._process_targets()

    def _process_targets(self):
        """Process and validate target specifications."""
        for target in self.targets:
            if 'weight' not in target:
                target['weight'] = 1.0
            if 'target_file' in target:
                import pandas as pd
                df = pd.read_csv(target['target_file'])
                if 'time' not in df.columns or 'value' not in df.columns:
                    raise ValueError(f"target_file must have 'time' and 'value' columns: {target['target_file']}")
                target['target_times'] = df['time'].values
                target['target_values'] = df['value'].values
            elif 'target_value' in target:
                target['target_value'] = float(target['target_value'])
            else:
                raise ValueError(f"Target must have either 'target_file' or 'target_value': {target}")

    def _get_simulated_value(
        self,
        target: Dict,
        simulated_values: Dict[str, Union[np.ndarray, float]]
    ) -> np.ndarray:
        """Look up and return simulated value for a target as numpy array."""
        name = target['name']
        extraction_type = target.get('type', 'time_series')
        target_key = f"{name}_{extraction_type}" if extraction_type != 'time_series' else name
        if target_key in simulated_values:
            sim_value = simulated_values[target_key]
        elif name in simulated_values:
            sim_value = simulated_values[name]
        else:
            raise ValueError(f"Simulated value for '{name}' (type: {extraction_type}) not found")
        if not isinstance(sim_value, np.ndarray):
            sim_value = np.array([sim_value]) if np.isscalar(sim_value) else np.array(sim_value)
        return np.asarray(sim_value)

    def _error_for_time_series(
        self,
        target: Dict,
        sim_value: np.ndarray,
        simulated_values: Dict
    ) -> float:
        """Compute error for a time series target."""
        name = target['name']
        target_times = np.array(target['target_times'])
        target_values = np.array(target['target_values'])
        sim_times = simulated_values.get(f'{name}_times')
        if sim_times is None:
            sim_times = np.linspace(0, 1, len(sim_value))
        sim_times = np.array(sim_times)

        if len(sim_times) != len(sim_value):
            raise ValueError(
                f"Time and value arrays have different lengths for {name}: "
                f"{len(sim_times)} vs {len(sim_value)}"
            )

        if len(sim_times) < 2 or len(target_times) == 0:
            return 1e10

        sim_interp = _interpolate_to_target_times(sim_times, sim_value, target_times)
        if not np.all(np.isfinite(sim_interp)):
            return 1e10

        return _l2_error(target_values, sim_interp, self.normalize)

    def _error_for_scalar(self, target: Dict, sim_value: np.ndarray) -> float:
        """Compute error for a scalar target (min, max, mean)."""
        target_value = float(target['target_value'])
        sim_scalar = float(sim_value.item() if sim_value.size == 1 else sim_value.flat[0])
        return _l2_error(
            np.array([target_value]), np.array([sim_scalar]), self.normalize
        )

    def compute(self, simulated_values: Dict[str, Union[np.ndarray, float]]) -> float:
        """
        Compute objective function value. Returns weighted sum of L2 error for all targets (both time series and scalars).

        Args:
            simulated_values: Dictionary mapping output names to simulated values

        Returns:
            Total weighted error
        """
        if self.custom_function:
            return self.custom_function(simulated_values, self.targets)

        total_error = 0.0
        for target in self.targets:
            sim_value = self._get_simulated_value(target, simulated_values)
            weight = float(target.get('weight', 1.0))
            extraction_type = target.get('type', 'time_series')

            if extraction_type == 'time_series':
                error = self._error_for_time_series(target, sim_value, simulated_values)
            else:
                error = self._error_for_scalar(target, sim_value)

            total_error += weight * error

        return float(total_error)


def create_objective(targets: List[Dict], **kwargs) -> ObjectiveFunction:
    """
    Create objective function. Use normalize=True for relative error, False for L2.
    """
    normalize = kwargs.pop('normalize', False)
    custom_function = kwargs.pop('custom_function', None)

    return ObjectiveFunction(
        targets=targets,
        normalize=normalize,
        custom_function=custom_function,
        **kwargs
    )
