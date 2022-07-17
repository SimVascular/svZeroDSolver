"""This module holds the BloodVesselJunction class."""
from .internaljunction import InternalJunction


class BloodVesselJunction(InternalJunction):
    """Blood vessel junction.

    (dummy for future implementation of blood pressure losses at junctions)

    Attributes:
        name: Name of the block.
        inflow_nodes: Inflow nodes of the element.
        outflow_nodes: Outflow nodes of the element.
    """

    pass
