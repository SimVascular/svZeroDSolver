from .block import Block
from .dofhandler import DOFHandler


class Node:
    """Node.

    Nodes connect two blocks with each other. Each node corresponds to a
    flow and pressure value of the system.

    Attributes:
        name: Name of the node.
        flow_dof: Global ID of the flow value associated with the node.
        pres_dof: Global ID of the pressure value associated with the node.
    """

    def __init__(self, ele1: Block, ele2: Block, name: str = None):
        """Create a new Node instance.

        Args:
            ele1: First element for the node to connect.
            ele2: Second element for the node to connect.
            name: Optional name of the node.
        """
        self.name = name
        self.flow_dof: int = None
        self.pres_dof: int = None

        # Make the node the ouflow node of the first element and the inflow
        # node of the second element
        ele1.outflow_nodes.append(self)
        ele2.inflow_nodes.append(self)

    def setup_dofs(self, dofhandler: DOFHandler):
        """Setup degree of freedoms of the node.

        Registers the pressure and the flow variable at a DOF handler.

        Args:
            dofhandler: The DOF handler to register the variables at.
        """
        self.flow_dof = dofhandler.register_variable("Q_" + self.name)
        self.pres_dof = dofhandler.register_variable("P_" + self.name)
