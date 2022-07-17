"""This module holds the InternalJunction class."""
import numpy as np

from .block import Block


class InternalJunction(Block):
    r"""Junction block.

    Represents a basic model junction without special mechanical behavior.
    Across all inlets and outlets of the junction, mass is conserved and
    pressure is continuous.

    Attributes:
        name: Name of the block.
        inflow_nodes: Inflow nodes of the element.
        outflow_nodes: Outflow nodes of the element.
    """

    def setup_dofs(self, dofhandler):
        """Setup degree of freedoms of the block.

        Registers equations and internal variables at a DOF handler.

        Args:
            dofhandler: The DOF handler to register the variables and equations
                at.
        """
        # Derive number of inlets and outlets
        num_inlets = len(self.inflow_nodes)
        num_outlets = len(self.outflow_nodes)

        # Set number of equations of a junction block based on number of
        # inlets/outlets. Must be set before calling parent constructor
        self._NUM_EQUATIONS = num_inlets + num_outlets
        super().setup_dofs(dofhandler)

        # Set some constant element element contributions that needed
        # _NUM_EQUATIONS
        self._mat["F"] = np.zeros(
            (self._NUM_EQUATIONS, self._NUM_EQUATIONS * 2)
        )
        for i in range(self._NUM_EQUATIONS - 1):
            self._mat["F"][i, [0, 2 * i + 2]] = [1.0, -1.0]
        self._mat["F"][-1, np.arange(1, 2 * self._NUM_EQUATIONS, 2)] = [
            1.0
        ] * num_inlets + [-1.0] * num_outlets

    @classmethod
    def from_config(cls, config: dict) -> "InternalJunction":
        """Create block from config dictionary.

        Args:
            config: The configuration dict for the block.

        Returns:
            The block instance.
        """
        name = config["junction_name"]
        if not name.startswith("J") and not name[1].isnumeric():
            raise ValueError(
                f"Invalid junction name {name}. Junction names must "
                "start with J following by a number."
            )
        return cls(name=name)
