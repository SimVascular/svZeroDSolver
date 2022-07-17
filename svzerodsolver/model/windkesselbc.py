"""This module holds the WindkesselBC class."""
from typing import Sequence

import numpy as np

from .block import Block


class WindkesselBC(Block):
    """Windkessel RCR boundary condition.

    Attributes:
        name: Name of the block.
        inflow_nodes: Inflow nodes of the element.
        outflow_nodes: Outflow nodes of the element.
    """

    _NUM_EQUATIONS = 2
    _NUM_INTERNAL_VARS = 1

    def __init__(self, params: dict = None, name: str = None):
        """Create a new WindkesselBC instance.

        Args:
            params: The configuration paramaters of the block. Mostly comprised
                of constants for element contribution calculation.
            name: Optional name of the block.
        """
        super().__init__(params=params, name=name)

        self._mat["E"] = np.zeros((2, 3), dtype=float)
        self._mat["F"] = np.array(
            [[1.0, 0.0, -1.0], [0.0, 0.0, -1.0]], dtype=float
        )
        self._vec["C"] = np.array([0.0, 0.0], dtype=float)

        self._rp_func = None
        if isinstance(self._params["Rp"], Sequence):
            self._rp_func = self._interpolate(
                self._params["time"], self._params["Rp"]
            )

        self._c_func = None
        if isinstance(self._params["C"], Sequence):
            self._c_func = self._interpolate(
                self._params["time"], self._params["C"]
            )

        self._rd_func = None
        if isinstance(self._params["Rd"], Sequence):
            self._rd_func = self._interpolate(
                self._params["time"], self._params["Rd"]
            )

        self._pd_func = None
        if isinstance(self._params["Pd"], Sequence):
            self._pd_func = self._interpolate(
                self._params["time"], self._params["Pd"]
            )

        if [self._rp_func, self._c_func, self._rd_func, self._pd_func] == [
            None
        ] * 4:
            self._mat["E"][1, 2] = -self._params["Rd"] * self._params["C"]
            self._mat["F"][0, 1] = -self._params["Rp"]
            self._mat["F"][1, 1] = self._params["Rd"]
            self._vec["C"][1] = self._params["Pd"]
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
            Rp=config["bc_values"].get("Rp"),
            C=config["bc_values"].get("C"),
            Rd=config["bc_values"].get("Rd"),
            Pd=config["bc_values"].get("Pd"),
        )
        return cls(params=params, name=config["name"])

    def update_time(self, time):
        """Update time dependent element contributions.

        Args:
            time: Current time.
        """
        rd_t = self._rd_func(time)
        self._mat["E"][1, 2] = -rd_t * self._c_func(time)
        self._mat["F"][0, 1] = -self._rp_func(time)
        self._mat["F"][1, 1] = rd_t
        self._vec["C"][1] = self._pd_func(time)
