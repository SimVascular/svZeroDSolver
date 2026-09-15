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
    User can provide:
    - single value
    - value + relative_bounds (percent)
    - target_range [min, max]
    Legacy alias: uncertainty
    """
    def _ordered_bounds(lo: np.ndarray, hi: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        lo_arr = np.asarray(lo, dtype=float)
        hi_arr = np.asarray(hi, dtype=float)
        return np.minimum(lo_arr, hi_arr), np.maximum(lo_arr, hi_arr)

    relative_bounds = target.get('relative_bounds', target.get('uncertainty'))
    if 'target_file' in target:
        # Time series
        t = np.asarray(target['target_values'])

        # A target_range is provided
        if 'target_range' in target:
            lo_val, hi_val = float(target['target_range'][0]), float(target['target_range'][1])
            return (np.full(len(t), lo_val), np.full(len(t), hi_val))

        # Relative bounds are provided as percent or [min, max]
        if relative_bounds is not None:
            pct = _parse_percent(relative_bounds)
            if pct is not None:
                return _ordered_bounds(t * (1.0 - pct), t * (1.0 + pct))
            if isinstance(relative_bounds, (list, tuple)) and len(relative_bounds) == 2:
                lo_val, hi_val = float(relative_bounds[0]), float(relative_bounds[1])
                return (np.full(len(t), lo_val), np.full(len(t), hi_val))
        
        # No target_range or relative_bounds, so point target
        return (t.copy(), t.copy())  # point target
    else:
        # Scalar

        # A target_range is provided
        if 'target_range' in target:
            lo_val, hi_val = float(target['target_range'][0]), float(target['target_range'][1])
            return (np.array([lo_val]), np.array([hi_val]))
        
        v = float(target['target_value'])
        
        # Relative bounds are provided as percent or [min, max]
        if relative_bounds is not None:
            pct = _parse_percent(relative_bounds)
            if pct is not None:
                return _ordered_bounds(
                    np.array([v * (1.0 - pct)]),
                    np.array([v * (1.0 + pct)]),
                )
            if isinstance(relative_bounds, (list, tuple)) and len(relative_bounds) == 2:
                return (np.array([float(relative_bounds[0])]), np.array([float(relative_bounds[1])]))
        
        # No target_range or relative_bounds, so point target
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


def _rel_errors_outside_range(
    sim_values: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
) -> np.ndarray:
    """
    Relative errors for values outside the target range [lo, hi].
    In-range values contribute zero; outside the range we use residual (distance to nearest bound).
    Works for both time series (length N) and scalar targets (length 1). For
    time series, each time point is treated as an individual scalar target.
    - if in range: residual = 0
    - if below lo: residual = lo - sim_value
    - if above hi: residual = sim_value - hi

    Residuals are then normalized by the midpoint.
    Returns: array of |residual_n| / |midpoint_n| per point, where midpoint_n = (lo_n + hi_n) / 2.
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
    return abs_residual / abs_midpoint


class ObjectiveFunction:
    """
    Objective function for optimization. Computes weighted error between
    simulated values and targets. Supports scalars and time series.

    Internally each target is stored as a range [lo, hi]. User can specify:
    - Single value (target_value or target_file): range = [v, v]
    - Value + relative_bounds percent: range = value * (1 ± pct)
    - target_range [min, max] directly
    Error is zero within the range; outside, penalty by distance from nearest bound.
    """

    def __init__(self, targets: List[Dict], norm: str):
        """
        Initialize objective function.

        Args:
            targets: List of target specifications
            norm: 'L1' for sum of absolute errors, 'L2' for Euclidean norm of error vector
        """
        if norm not in ("L1", "L2"):
            raise ValueError(
                f"norm is required and must be 'L1' or 'L2', got {norm!r}. "
                "Options: L1 = sum of absolute relative errors; L2 = Euclidean norm of the error vector. "
                "Specify in your tuning YAML under objective: { norm: L1 } or { norm: L2 }."
            )
        self.targets = targets
        self.norm = norm
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
    ) -> np.ndarray:
        """Return array of relative errors, one per time point."""
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
            return np.array([1e10])

        sim_interp = _interpolate_to_target_times(sim_times, sim_value, target_times)
        if not np.all(np.isfinite(sim_interp)):
            return np.array([1e10])

        return _rel_errors_outside_range(
            sim_interp, target['range_lo'], target['range_hi']
        )

    def _errors_for_scalar(self, target: Dict, sim_value: np.ndarray) -> np.ndarray:
        """Return array of one relative error for a scalar target."""
        sim_scalar = float(sim_value.item() if sim_value.size == 1 else sim_value.flat[0])
        return _rel_errors_outside_range(
            np.array([sim_scalar]), target['range_lo'], target['range_hi']
        )

    def compute(self, simulated_values: Dict[str, Union[np.ndarray, float]]) -> float:
        """
        Compute objective function value. Collects weighted relative errors from all
        targets (each time series point or scalar contributes one error value), then returns
        L1 norm (sum of absolute errors) or L2 norm (Euclidean norm) of that vector.

        Args:
            simulated_values: Dictionary mapping output names to simulated values

        Returns:
            Total error: L1 or L2 norm of the weighted relative-error vector
        """
        error_chunks: List[np.ndarray] = []
        for target in self.targets:
            sim_value = self._get_simulated_value(target, simulated_values)
            weight = float(target['weight'])
            target_type = target['type']

            if target_type == 'time_series':
                errors = self._error_for_time_series(target, sim_value, simulated_values)
            else:
                errors = self._errors_for_scalar(target, sim_value)

            error_chunks.append(weight * errors)
        all_errors = np.concatenate(error_chunks) if error_chunks else np.array([])

        if self.norm == "L1":
            ord = 1
        elif self.norm == "L2":
            ord = 2
        else:
            raise ValueError(f"norm must be 'L1' or 'L2', got {self.norm!r}")
        return float(np.linalg.norm(all_errors, ord=ord))


def create_objective(targets: List[Dict], **kwargs) -> ObjectiveFunction:
    """
    Create objective function object.
    Targets with 'relative_bounds' (or legacy 'uncertainty') or target_range use
    range-based penalty.
    """
    norm = kwargs.pop("norm")
    return ObjectiveFunction(targets=targets, norm=norm, **kwargs)
