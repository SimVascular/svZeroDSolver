from typing import Sequence

import numpy as np

from .block import Block


class OpenLoopCoronaryBC(Block):
    """
    open-loop coronary BC = RCRCR BC
    Publication reference: Kim, H. J. et al. Patient-specific modeling of blood flow and pressure in human coronary arteries. Annals of Biomedical Engineering 38, 3195–3209 (2010)."
    """

    _NUM_EQUATIONS = 2
    _NUM_INTERNAL_VARS = 1

    def __init__(self, params: dict = None, name: str = None):
        super().__init__(params=params, name=name)

        self._pv_func = None
        if isinstance(self._params["Pv"], Sequence):
            self._pv_func = self._interpolate(self._params["time"], self._params["Pv"])

        self._pim_func = None
        if isinstance(self._params["Pim"], Sequence):
            self._pim_func = self._interpolate(
                self._params["time"], self._params["Pim"]
            )

        self._vec["C"] = np.zeros(2)
        self._mat["F"] = np.zeros((2, 3))

        self._mat["F"][0, 2] = -1.0

        Cim_Rv = self._params["Cim"] * self._params["Rv"]
        self._mat["E"] = np.zeros((2, 3))
        self._mat["E"][0, 0] = -self._params["Ca"] * Cim_Rv
        self._mat["E"][0, 1] = self._params["Ra"] * self._params["Ca"] * Cim_Rv
        self._mat["E"][0, 2] = -Cim_Rv
        self._mat["E"][1, 2] = -Cim_Rv * self._params["Ram"]
        self._mat["F"][0, 1] = Cim_Rv
        self._mat["F"][1, 0] = Cim_Rv
        self._mat["F"][1, 1] = -Cim_Rv * self._params["Ra"]
        self._mat["F"][1, 2] = -(self._params["Rv"] + self._params["Ram"])

        if self._pv_func is None and self._pim_func is None:
            self._vec["C"][0] = (
                -self._params["Cim"] * self._params["Pim"]
                + self._params["Cim"] * self._params["Pv"]
            )
            # Pa is assumed to be 0.0
            self._vec["C"][1] = (
                -self._params["Cim"]
                * (self._params["Rv"] + self._params["Ram"])
                * self._params["Pim"]
                + self._params["Ram"] * self._params["Cim"] * self._params["Pv"]
            )
            self.update_time = super().update_time
        elif self._pv_func is None:
            self._pv_func = lambda _: self._params["Pv"]
        elif self._pim_func is None:
            self._pim_func = lambda _: self._params["Pim"]

        if self._params["steady"]:
            del self._mat["E"]
            self._mat["F"] = np.array(
                [
                    [
                        -self._params["Cim"],
                        self._params["Cim"]
                        * (self._params["Ra"] + self._params["Ram"]),
                        1.0,
                    ],
                    [
                        -1.0,
                        self._params["Ra"] + self._params["Ram"] + self._params["Rv"],
                        0.0,
                    ],
                ]
            )
            self._vec["C"][0] = -self._params["Cim"] * self._params["Pim"]
            self._vec["C"][1] = self._params["Pv"]

    @classmethod
    def from_config(cls, config):
        params = dict(
            time=config["bc_values"].get("t", None),
            Ra=config["bc_values"]["Ra1"],
            Ca=config["bc_values"]["Ca"],
            Ram=config["bc_values"]["Ra2"],
            Cim=config["bc_values"]["Cc"],
            Rv=config["bc_values"]["Rv1"],
            Pim=config["bc_values"]["Pim"],
            Pv=config["bc_values"]["P_v"],
            steady=config["steady"],
        )

        return cls(params=params, name=config["name"])

    def update_time(self, time):
        # For this open-loop coronary BC, the ordering of solution unknowns is : (P_in, Q_in, V_im)
        # where V_im is the volume of the second capacitor, Cim
        # Q_in is the flow through the first resistor
        # and P_in is the pressure at the inlet of the first resistor
        Pim_value = self._pim_func(time)
        Pv_value = self._pv_func(time)
        self._vec["C"][0] = (
            -self._params["Cim"] * Pim_value + self._params["Cim"] * Pv_value
        )
        # Pa is assumed to be 0.0
        self._vec["C"][1] = (
            -self._params["Cim"]
            * (self._params["Rv"] + self._params["Ram"])
            * Pim_value
            + self._params["Ram"] * self._params["Cim"] * Pv_value
        )
