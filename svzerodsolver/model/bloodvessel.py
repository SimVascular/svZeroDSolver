# Copyright (c) Stanford University, The Regents of the University of
#               California, and others.
#
# All Rights Reserved.
#
# See Copyright-SimVascular.txt for additional details.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""This module holds the BloodVessel class."""
import numpy as np

from .block import Block


class BloodVessel(Block):
    r"""Resistor-capacitor-inductor blood vessel with optional stenosis.

    Models the mechanical behavior of a bloodvessel with optional stenosis.

    Attributes:
        name: Name of the block.
        inflow_nodes: Inflow nodes of the element.
        outflow_nodes: Outflow nodes of the element.
    """

    _NUM_EQUATIONS = 3
    _NUM_INTERNAL_VARS = 1

    def __init__(self, params: dict = None, name: str = None):
        """Create a new BloodVessel instance.

        Args:
            params: The configuration paramaters of the block. Mostly comprised
                of constants for element contribution calculation.
            name: Optional name of the block.
        """
        super().__init__(params=params, name=name)

        # the ordering of the solution variables is : (P_in, Q_in, P_out, Q_out)
        self._mat["E"] = np.zeros((3, 5), dtype=float)
        self._mat["E"][0, 3] = -self._params["L"]
        self._mat["E"][1, 4] = -self._params["C"]
        self._mat["F"] = np.array(
            [
                [1.0, -self._params["R"], -1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, -1.0, 0.0],
                [1.0, -self._params["R"], 0.0, 0.0, -1.0],
            ],
            dtype=float,
        )
        self._mat["dF"] = np.zeros((3, 5), dtype=float)

        # This class's update_solution is only needed if class is configured with
        # non-zero stenosis coefficient. If stenosis coefficient is 0.0 the
        # emtpy parent update_solution can be used and is faster.
        if self._params["stenosis_coefficient"] == 0.0:
            self.update_solution = super().update_solution

    @classmethod
    def from_config(cls, config: dict) -> "BloodVessel":
        """Create block from config dictionary.

        Args:
            config: The configuration dict for the block.

        Returns:
            The block instance.
        """
        params = dict(
            R=config["zero_d_element_values"].get("R_poiseuille"),
            C=config["zero_d_element_values"].get("C", 0.0),
            L=config["zero_d_element_values"].get("L", 0.0),
            stenosis_coefficient=config["zero_d_element_values"].get(
                "stenosis_coefficient", 0.0
            ),
        )
        return cls(params=params, name="V" + str(config["vessel_id"]))

    def update_solution(self, y: np.ndarray) -> None:
        """Update solution dependent element contributions.

        Args:
            y: Current solution.
        """
        Q_in = np.abs(y[self.inflow_nodes[0].flow_dof])
        fac1 = -self._params["stenosis_coefficient"] * Q_in
        fac2 = fac1 - self._params["R"]
        self._mat["F"][[0, 2], 1] = fac2
        self._mat["dF"][[0, 2], 1] = fac1
