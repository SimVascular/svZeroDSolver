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
"""This module holds the Node class."""
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
