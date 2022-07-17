"""This module holds the ResistanceBC class."""
from typing import Sequence

import numpy as np

from .block import Block


class ResistanceBC(Block):
    """Resistance boundary condition.

    Attributes:
        name: Name of the block.
        inflow_nodes: Inflow nodes of the element.
        outflow_nodes: Outflow nodes of the element.
    """

    _NUM_EQUATIONS = 1

    def __init__(self, params: dict = None, name: str = None):
        """Create a new ResistanceBC instance.

        Args:
            params: The configuration paramaters of the block. Mostly comprised
                of constants for element contribution calculation.
            name: Optional name of the block.
        """
        super().__init__(params=params, name=name)

        self._r_func = None
        if isinstance(self._params["R"], Sequence):
            self._r_func = self._interpolate(
                self._params["time"], self._params["R"]
            )
            self._mat["F"] = np.array([[1.0, 0.0]], dtype=float)
        else:
            self._mat["F"] = np.array([[1.0, -self._params["R"]]], dtype=float)

        self._pd_func = None
        if isinstance(self._params["Pd"], Sequence):
            self._pd_func = self._interpolate(
                self._params["time"], self._params["Pd"]
            )
            self._vec["C"] = np.array([0.0], dtype=float)
        else:
            self._vec["C"] = np.array([-self._params["Pd"]], dtype=float)

        if self._r_func is None and self._pd_func is None:
            self.update_time = super().update_time

    @classmethod
    def from_config(cls, config):
        """Create block from config dictionary.

        Args:
            config: The configuration dict for the block.

        Returns:
            The block instance.
        """
        params = dict(
            time=config["bc_values"].get("t", None),
            Pd=config["bc_values"].get("Pd"),
            R=config["bc_values"].get("R"),
        )
        return cls(params=params, name=config["name"])

    def update_time(self, time):
        """Update time dependent element contributions.

        Args:
            time: Current time.
        """
        self._mat["F"][0, 1] = -self._r_func(time)
        self._vec["C"][0] = -self._pd_func(time)
