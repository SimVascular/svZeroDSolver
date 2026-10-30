"""
Parameter scaling for optimization.

Scaling objects convert between physical (model) space and optimizer space.
E.g. log scaling: optimizer works with log(p); to_physical applies exp, to_opt applies log.
"""

import numpy as np
from typing import List


class ParameterScaling:
    """
    Converts a single parameter between physical space and optimizer space.
    """

    @property
    def requires_positive(self) -> bool:
        """If True, physical values and bounds must be strictly positive."""
        return False

    def to_physical(self, x_opt: float) -> float:
        """Convert one value from optimizer space to physical space."""
        raise NotImplementedError

    def to_opt(self, x_physical: float) -> float:
        """Convert one value from physical space to optimizer space."""
        raise NotImplementedError


class IdentityScaling(ParameterScaling):
    """No transformation; physical and optimizer space are the same."""

    def to_physical(self, x_opt: float) -> float:
        return float(x_opt)

    def to_opt(self, x_physical: float) -> float:
        return float(x_physical)


class LogScaling(ParameterScaling):
    """Optimizer works with log(parameter); physical = exp(opt)."""

    @property
    def requires_positive(self) -> bool:
        return True

    def to_physical(self, x_opt: float) -> float:
        return float(np.exp(x_opt))

    def to_opt(self, x_physical: float) -> float:
        return float(np.log(x_physical))


class MaxScaling(ParameterScaling):
    """Scale by the max absolute bound magnitude: optimizer space = physical / max(|lo|, |hi|)."""

    def __init__(self, bounds: tuple):
        lo, hi = float(bounds[0]), float(bounds[1])
        self._max_bound = max(abs(lo), abs(hi))
        if self._max_bound == 0:
            raise ValueError("Max scaling requires non-zero bounds")

    def to_physical(self, x_opt: float) -> float:
        return float(x_opt * self._max_bound)

    def to_opt(self, x_physical: float) -> float:
        return float(x_physical / self._max_bound)


def get_scaling(name: str, bounds: tuple = None) -> ParameterScaling:
    """
    Return a scaling instance from a config string.

    Args:
        name: Scaling name: 'identity', 'log', or 'max'. 'identity' or None = no transform.
        bounds: Optional (lo, hi) for the parameter. Required when name is 'max'.

    Returns:
        A ParameterScaling instance.
    """
    if name == "identity" or name is None:
        return IdentityScaling()
    if name == "log":
        return LogScaling()
    if name == "max":
        if bounds is None or len(bounds) != 2:
            raise ValueError("Scaling 'max' requires bounds [lo, hi]")
        return MaxScaling(tuple(bounds))
    raise ValueError(f"Unknown scaling {name!r}; use 'identity', 'log', or 'max'")


def to_physical_array(x_opt: np.ndarray, scalings: List[ParameterScaling]) -> np.ndarray:
    """
    Convert a full parameter vector from optimizer space to physical space.

    Args:
        x_opt: Parameter vector in optimizer space
        scalings: List of ParameterScaling instances for each parameter

    Returns:
        Parameter vector in physical space
    """
    out = np.empty_like(x_opt, dtype=float)
    for i, s in enumerate(scalings):
        out[i] = s.to_physical(x_opt[i])
    return out


def to_opt_array(x_physical: np.ndarray, scalings: List[ParameterScaling]) -> np.ndarray:
    """
    Convert a full parameter vector from physical space to optimizer space.

    Args:
        x_physical: Parameter vector in physical space
        scalings: List of ParameterScaling instances for each parameter

    Returns:
        Parameter vector in optimizer space
    """
    out = np.empty_like(x_physical, dtype=float)
    for i, s in enumerate(scalings):
        out[i] = s.to_opt(x_physical[i])
    return out
