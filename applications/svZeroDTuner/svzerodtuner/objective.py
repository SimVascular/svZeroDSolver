"""
Objective function for sv0D Tuning Framework.
"""

import numpy as np
from typing import Dict, List, Optional, Union, Tuple
from scipy.interpolate import interp1d


def _parse_percent(x) -> Optional[float]:
    """Parse percent: 5, '5%' -> 0.05. Returns None if not a percent form."""
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x) / 100.0
    if isinstance(x, str) and x.strip().endswith('%'):
        try:
            return float(x.strip()[:-1]) / 100.0
        except ValueError:
            return None
    return None


def _compute_range(target: Dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert user spec to (lo, hi) range. Internally we only keep range.
    User can provide: single value, value+uncertainty (percent), or target_range [min,max].
    """
    if 'target_file' in target:
        # Time series
        t = np.asarray(target['target_values'])

        # A target_range is provided
        if 'target_range' in target:
            lo_val, hi_val = float(target['target_range'][0]), float(target['target_range'][1])
            return (np.full(len(t), lo_val), np.full(len(t), hi_val))

        # An uncertainty percentage is provided
        if 'uncertainty' in target:
            unc = target['uncertainty']
            pct = _parse_percent(unc)
            if pct is not None:
                return (t * (1.0 - pct), t * (1.0 + pct))
            if isinstance(unc, (list, tuple)) and len(unc) == 2:
                lo_val, hi_val = float(unc[0]), float(unc[1])
                return (np.full(len(t), lo_val), np.full(len(t), hi_val))
        
        # No target_range or uncertainty, so point target
        return (t.copy(), t.copy())  # point target
    else:
        # Scalar

        # A target_range is provided
        if 'target_range' in target:
            lo_val, hi_val = float(target['target_range'][0]), float(target['target_range'][1])
            return (np.array([lo_val]), np.array([hi_val]))
        
        v = float(target['target_value'])
        
        # An uncertainty percentage is provided
        if 'uncertainty' in target:
            unc = target['uncertainty']
            pct = _parse_percent(unc)
            if pct is not None:
                return (np.array([v * (1.0 - pct)]), np.array([v * (1.0 + pct)]))
            if isinstance(unc, (list, tuple)) and len(unc) == 2:
                return (np.array([float(unc[0])]), np.array([float(unc[1])]))
        
        # No target_range or uncertainty, so point target
        return (np.array([v]), np.array([v]))  


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


def _sum_rel_errors_outside_range(
    sim_values: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
) -> float:
    """
    Sum of relative errors for values outside the target range [lo, hi].
    In-range values contribute zero; outside the range we use residual (distance to nearest bound).
    Works for both time series (length N) and scalar targets (length 1). For
    time series, each time point is treated as an individual scalar target.
    - if in range: residual = 0
    - if below lo: residual = lo - sim_value
    - if above hi: residual = sim_value - hi

    Residuals are then normalized by the midpoint, per time point.

    Returns: Σ |residual_n| / |midpoint_n|, where midpoint_n = (lo_n + hi_n) / 2. 
    """
    sim_values = np.asarray(sim_values)
    lo = np.asarray(lo)
    hi = np.asarray(hi)
    below = sim_values < lo
    above = sim_values > hi
    residual = np.zeros_like(sim_values, dtype=float)
    residual[below] = lo[below] - sim_values[below]
    residual[above] = sim_values[above] - hi[above]
    abs_residual = np.abs(residual)
    midpoint = (lo + hi) / 2.0
    abs_midpoint = np.maximum(np.abs(midpoint), 1e-14)
    return float(np.sum(abs_residual / abs_midpoint))


class ObjectiveFunction:
    """
    Objective function for optimization. Computes weighted error between
    simulated values and targets. Supports scalars and time series.

    Internally each target is stored as a range [lo, hi]. User can specify:
    - Single value (target_value or target_file): range = [v, v]
    - Value + uncertainty percent: range = value * (1 ± pct)
    - target_range [min, max] directly
    Error is zero within the range; outside, penalty by distance from nearest bound.
    """

    def __init__(self, targets: List[Dict]):
        """
        Initialize objective function.

        Args:
            targets: List of target specifications
        """
        self.targets = targets
        self._process_targets()

    def _process_targets(self):
        """Process targets: load data and compute range [lo, hi] for each."""
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
            elif 'target_range' not in target:
                raise ValueError(f"Target must have 'target_file', 'target_value', or 'target_range': {target}")
            lo, hi = _compute_range(target)
            target['range_lo'] = lo
            target['range_hi'] = hi

    def _get_simulated_value(
        self,
        target: Dict,
        simulated_values: Dict[str, Union[np.ndarray, float]]
    ) -> np.ndarray:
        """Look up and return simulated value for a target as numpy array."""
        name = target['name']
        if name not in simulated_values:
            raise ValueError(f"Simulated value for '{name}' not found")
        sim_value = simulated_values[name]
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

        return _sum_rel_errors_outside_range(
            sim_interp, target['range_lo'], target['range_hi']
        )

    def _error_for_scalar(self, target: Dict, sim_value: np.ndarray) -> float:
        """Compute error for a scalar target."""
        sim_scalar = float(sim_value.item() if sim_value.size == 1 else sim_value.flat[0])
        return _sum_rel_errors_outside_range(
            np.array([sim_scalar]), target['range_lo'], target['range_hi']
        )

    def compute(self, simulated_values: Dict[str, Union[np.ndarray, float]]) -> float:
        """
        Compute objective function value. Returns weighted sum of relative errors 
        for all targets. For time series, each time point is treated as an 
        individual scalar target.

        Args:
            simulated_values: Dictionary mapping output names to simulated values

        Returns:
            Total weighted error
        """
        total_error = 0.0
        for target in self.targets:
            sim_value = self._get_simulated_value(target, simulated_values)
            weight = float(target.get('weight', 1.0))
            target_type = target.get('type', 'time_series')

            if target_type == 'time_series':
                error = self._error_for_time_series(target, sim_value, simulated_values)
            else:
                error = self._error_for_scalar(target, sim_value)

            total_error += weight * error

        return float(total_error)


def create_objective(targets: List[Dict], **kwargs) -> ObjectiveFunction:
    """
    Create objective function object.
    Targets with 'uncertainty' (percent or [min,max]) or target_range use range-based penalty.
    """
    return ObjectiveFunction(targets=targets, **kwargs)
