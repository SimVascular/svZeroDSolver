from typing import Sequence

import numpy as np

from .block import Block


class PressureReferenceBC(Block):
    """
    Unsteady P reference
    """

    _NUM_EQUATIONS = 1

    def __init__(self, params: dict = None, name: str = None):
        super().__init__(params=params, name=name)

        if isinstance(self._params["P"], Sequence):
            self._p_func = self._interpolate(self._params["time"], self._params["P"])
            self._vec["C"] = np.zeros(1, dtype=float)
        else:
            self._vec["C"] = np.array([-self._params["P"]], dtype=float)
            self.update_time = super().update_time

        self._mat["F"] = np.array([[1.0, 0.0]], dtype=float)

    @classmethod
    def from_config(cls, config):
        params = dict(
            time=config["bc_values"].get("t", None),
            P=config["bc_values"].get("P"),
        )
        return cls(params=params, name=config["name"])

    def update_time(self, time):
        self._vec["C"][0] = -self._p_func(time)
