from typing import Sequence

import numpy as np

from .block import Block


class FlowReferenceBC(Block):
    """
    Flow reference
    """

    _NUM_EQUATIONS = 1

    def __init__(self, params: dict = None, name: str = None):
        super().__init__(params=params, name=name)

        if isinstance(self._params["Q"], Sequence):
            self._q_func = self._interpolate(self._params["time"], self._params["Q"])
            self._vec["C"] = np.zeros(1, dtype=float)
        else:
            self._vec["C"] = np.array([-self._params["Q"]], dtype=float)
            self.update_time = super().update_time

        self._mat["F"] = np.array([[0.0, 1.0]], dtype=float)

    @classmethod
    def from_config(cls, config):
        params = dict(
            time=config["bc_values"].get("t", None),
            Q=config["bc_values"].get("Q"),
        )
        return cls(params=params, name=config["name"])

    def update_time(self, time):
        self._vec["C"][0] = -self._q_func(time)
