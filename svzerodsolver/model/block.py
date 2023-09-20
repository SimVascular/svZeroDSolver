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
"""This module holds the Block class."""
from __future__ import annotations
from abc import ABC, abstractclassmethod

import numpy as np
from scipy.interpolate import CubicSpline

from .dofhandler import DOFHandler


class Block(ABC):
    """Base class for 0D model components.

    A block is the base class of 0D model elements. It is the place where the
    contribution of an element to the global system is controlled. It defines
    all relevant attributes and methods of an element and a few common helpers
    for setting it up.

    Attributes:
        name: Name of the block.
        inflow_nodes: Inflow nodes of the element.
        outflow_nodes: Outflow nodes of the element.
    """

    # Number of equations the element contirbutes at assembly
    _NUM_EQUATIONS = None

    # Number of internal variables of the block
    _NUM_INTERNAL_VARS = 0

    def __init__(self, params: dict = None, name: str = None):
        """Create a new Block.

        Args:
            params: The configuration paramaters of the block. Mostly comprised
                of constants for element contribution calculation.
            name: Optional name of the block.
        """
        self.name = name
        self._params = params

        # Inflow and outflow wires of the block. Will be set by wires
        self.inflow_nodes = []
        self.outflow_nodes = []

        # Element contribution matrices
        self._mat = {}
        self._vec = {}

        # row and column indices of block in global matrix
        self._global_row_id = None
        self._flat_row_ids = None
        self._flat_col_ids = None

    @abstractclassmethod
    def from_config(cls, config: dict) -> "Block":
        """Create block from config dictionary.

        Args:
            config: The configuration dict for the block.

        Returns:
            The block instance.
        """
        pass

    def setup_dofs(self, dofhandler: DOFHandler) -> None:
        """Setup degree of freedoms of the block.

        Registers equations and internal variables at a DOF handler.

        Args:
            dofhandler: The DOF handler to register the variables and equations
                at.
        """

        # Register internal variables
        internal_vars = [
            dofhandler.register_variable(f"var_{i}_{self.name}")
            for i in range(self._NUM_INTERNAL_VARS)
        ]

        # Collect assembly column ids based from inflow/outflow wires and
        # internal variables
        global_col_id = []
        for wire in self.inflow_nodes:
            global_col_id += [wire.pres_dof, wire.flow_dof]
        for wire in self.outflow_nodes:
            global_col_id += [wire.pres_dof, wire.flow_dof]
        global_col_id += internal_vars

        # Register equations of the block
        self._global_row_id = [
            dofhandler.register_equation() for _ in range(self._NUM_EQUATIONS)
        ]

        # Create flat indices to assemble matrices as flattend array (faster)
        meshgrid = np.array(
            np.meshgrid(self._global_row_id, global_col_id)
        ).T.reshape(-1, 2)
        self._flat_row_ids, self._flat_col_ids = meshgrid[:, 0], meshgrid[:, 1]

    def assemble(self, mat: dict[str, np.ndarray]) -> None:
        """Assemble block to global system.

        Args:
            mat: Global system.
        """
        for key, value in self._vec.items():
            mat[key][self._global_row_id] = value
        for key, value in self._mat.items():
            mat[key][self._flat_row_ids, self._flat_col_ids] = value.ravel()

    def update_time(self, time: float) -> None:
        """Update time dependent element contributions.

        Args:
            time: Current time.
        """
        pass

    def update_solution(self, y: np.ndarray) -> None:
        """Update solution dependent element contributions.

        Args:
            y: Current solution.
        """
        pass

    def _interpolate(self, times, values):
        if times is None:
            raise ValueError("No time sequence provided for interpolation.")
        return CubicSpline(
            np.array(times), np.array(values), bc_type="periodic"
        )
