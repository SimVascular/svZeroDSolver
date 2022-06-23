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

from abc import ABC, abstractclassmethod
import numpy as np

from scipy.interpolate import CubicSpline

from typing import Sequence


class DOFHandler:
    """Degree-of-freedom handler.

    Class for handling global degrees of freedom of a system. Assigns global
    IDs to variables and equations.

    Attributes:
        variables: List of variable names corresponding to the global IDs.
            Variables without a name have an entry None.
        n: Size of the system.
    """

    def __init__(self):
        """Create a new DOFHandler instance."""
        self._var_counter = -1
        self._eqn_counter = -1
        self.variables: list[str] = []

    @property
    def n(self) -> int:
        """Size of the system."""
        return self._eqn_counter + 1

    def register_variable(self, name: str = None) -> int:
        """Register a new variable and get global index.

        Args:
            name: Optional name of the variables.

        Returns:
            global_id: Global id of the variable.
        """
        self._var_counter += 1
        self.variables.append(name)
        return self._var_counter

    def register_equation(self) -> int:
        """Register a new equation and get global index.

        Returns:
            global_id: Global id of the equation.
        """
        self._eqn_counter += 1
        return self._eqn_counter


class Block(ABC):
    """Base class for lumped-parameter blocks.

    A block stores all information about the mechanical characteristics of
    a lumped-parameter element. The block is used to setup, update and
    assemble element contributions in a network.

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
        self.inflow_nodes: list[Node] = []
        self.outflow_nodes: list[Node] = []

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
        meshgrid = np.array(np.meshgrid(self._global_row_id, global_col_id)).T.reshape(
            -1, 2
        )
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

    def _interpolate(self, times, values, method="cubic_spline"):
        if times is None:
            raise ValueError("No time sequence provided for interpolation.")
        return CubicSpline(np.array(times), np.array(values), bc_type="periodic")


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


class Junction(Block):
    r"""Junction block.

    Represents a basic model junction without special mechanical behavior.
    Across all inlets and outlets of the junction, mass is conserved and
    pressure is continuous. The block contributes as many equations to the
    global system as it has inlet/outlet nodes :math:`n_{eq}=n_{inlets}+n_{outlets}`.
    One equation comes from the mass conservation:

    .. math::

        \sum_{i}^{n_{inlets}} Q_{i}=\sum_{j}^{n_{outlets}} Q_{j}

    And the remaining :math:`n_{eq}-1` equations come from the continuous
    pressure condition between each pair of different pressure values.

    .. math::

        P_{i}=P_{j} \quad \text{with} \quad i \neq j

    The ordering of the variables is:

    .. math::

        \left[\begin{array}{llllllllll}P_{\text {in}, 1}^{e} & Q_{\text {in}, 1}^{e} & \dots & P_{\text {in}, i}^{e} & Q_{\text {in}, i}^{e} & P_{\text {out }, 1}^{e} & Q_{\text {out }, 1}^{e} & \dots & P_{\text {out}, i}^{e} & Q_{\text {out}, i}^{e}\end{array}\right]
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
        self._mat["F"] = np.zeros((self._NUM_EQUATIONS, self._NUM_EQUATIONS * 2))
        for i in range(self._NUM_EQUATIONS - 1):
            self._mat["F"][i, [0, 2 * i + 2]] = [1.0, -1.0]
        self._mat["F"][-1, np.arange(1, 2 * self._NUM_EQUATIONS, 2)] = [
            1.0
        ] * num_inlets + [-1.0] * num_outlets

    @classmethod
    def from_config(cls, config: dict) -> "Junction":
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


class BloodVessel(Block):
    r"""Resistor-capacitor-inductor blood vessel with optional stenosis.

    Valid parameters:
        * :code:`R`: Poiseuille resistance :math:`\frac{8 \mu L}{\pi r^{4}}`
        * :code:`L`: Inductance
        * :code:`C`: Capacitance
        * :code:`stenosis_coefficient`: Stenosis coefficient :math:`K_{t} \frac{\rho}{2 A_{o}^{2}}\left(\frac{A_{o}}{A_{s}}-1\right)^{2}`.

    The governing equations for the local resistor-capacitor-inductor element are

    .. math::

        P_{i n}^{e}-P_{o u t}^{e}-R Q_{i n}^{e}-L \frac{d Q_{o u t}^{e}}{d t}=0

    .. math::

        Q_{\text {in }}^{e}-Q_{\text {out }}^{e}-C \frac{d P_{c}^{e}}{d t}=0

    .. math::

        P_{i n}^{e}-R Q_{i n}^{e}-P_{c}=0

    Stenosis:
        equation: delta_P = ( K_t * rho / ( 2 * (A_0)**2 ) ) * ( ( A_0 / A_s ) - 1 )**2 * Q * abs(Q) + R_poiseuille * Q
                          =               stenosis_coefficient                          * Q * abs(Q) + R_poiseuille * Q

        source: Mirramezani, M., Shadden, S.C. A distributed lumped parameter model of blood flow. Annals of Biomedical Engineering. 2020.
    """

    _NUM_EQUATIONS = 3
    _NUM_INTERNAL_VARS = 1

    def __init__(self, params: dict = None, name: str = None):
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


class ResistanceWithDistalPressure(Block):

    _NUM_EQUATIONS = 1

    def __init__(self, params: dict = None, name: str = None):
        super().__init__(params=params, name=name)

        self._r_func = None
        if isinstance(self._params["R"], Sequence):
            self._r_func = self._interpolate(self._params["time"], self._params["R"])
            self._mat["F"] = np.array([[1.0, 0.0]], dtype=float)
        else:
            self._mat["F"] = np.array([[1.0, -self._params["R"]]], dtype=float)

        self._pd_func = None
        if isinstance(self._params["Pd"], Sequence):
            self._pd_func = self._interpolate(self._params["time"], self._params["Pd"])
            self._vec["C"] = np.array([0.0], dtype=float)
        else:
            self._vec["C"] = np.array([-self._params["Pd"]], dtype=float)

        if self._r_func is None and self._pd_func is None:
            self.update_time = super().update_time

    @classmethod
    def from_config(cls, config):
        params = dict(
            time=config["bc_values"].get("t", None),
            Pd=config["bc_values"].get("Pd"),
            R=config["bc_values"].get("R"),
        )
        return cls(params=params, name=config["name"])

    def update_time(self, time):
        """
        the ordering is : (P_in,Q_in)
        """
        self._mat["F"][0, 1] = -self._r_func(time)
        self._vec["C"][0] = -self._pd_func(time)


class PressureRef(Block):
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


class FlowRef(Block):
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


class RCRBlockWithDistalPressure(Block):
    """
    Unsteady RCR - time-varying RCR values
    Formulation includes additional variable : internal pressure proximal to capacitance.
    """

    _NUM_EQUATIONS = 2
    _NUM_INTERNAL_VARS = 1

    def __init__(self, params: dict = None, name: str = None):
        super().__init__(params=params, name=name)

        self._mat["E"] = np.zeros((2, 3), dtype=float)
        self._mat["F"] = np.array([[1.0, 0.0, -1.0], [0.0, 0.0, -1.0]], dtype=float)
        self._vec["C"] = np.array([0.0, 0.0], dtype=float)

        self._rp_func = None
        if isinstance(self._params["Rp"], Sequence):
            self._rp_func = self._interpolate(self._params["time"], self._params["Rp"])

        self._c_func = None
        if isinstance(self._params["C"], Sequence):
            self._c_func = self._interpolate(self._params["time"], self._params["C"])

        self._rd_func = None
        if isinstance(self._params["Rd"], Sequence):
            self._rd_func = self._interpolate(self._params["time"], self._params["Rd"])

        self._pd_func = None
        if isinstance(self._params["Pd"], Sequence):
            self._pd_func = self._interpolate(self._params["time"], self._params["Pd"])

        if [self._rp_func, self._c_func, self._rd_func, self._pd_func] == [None] * 4:
            self._mat["E"][1, 2] = -self._params["Rd"] * self._params["C"]
            self._mat["F"][0, 1] = -self._params["Rp"]
            self._mat["F"][1, 1] = self._params["Rd"]
            self._vec["C"][1] = self._params["Pd"]
            self.update_time = super().update_time

    @classmethod
    def from_config(cls, config):
        params = dict(
            time=config["bc_values"].get("t", None),
            Rp=config["bc_values"].get("Rp"),
            C=config["bc_values"].get("C"),
            Rd=config["bc_values"].get("Rd"),
            Pd=config["bc_values"].get("Pd"),
        )
        return cls(params=params, name=config["name"])

    def update_time(self, time):
        """
        unknowns = [P_in, Q_in, internal_var (Pressure at the intersection of the Rp, Rd, and C elements)]
        """
        rd_t = self._rd_func(time)
        self._mat["E"][1, 2] = -rd_t * self._c_func(time)
        self._mat["F"][0, 1] = -self._rp_func(time)
        self._mat["F"][1, 1] = rd_t
        self._vec["C"][1] = self._pd_func(time)


class OpenLoopCoronaryWithDistalPressureBlock(Block):
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


def create_junction_blocks(parameters):
    """
    Purpose:
        Create the junction blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
    Returns:
        void, but updates parameters["blocks"] to include the junction_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    block_dict = {}
    connections = []
    for config in parameters["junctions"]:
        if config["junction_type"] in ["NORMAL_JUNCTION", "internal_junction"]:
            junction = Junction.from_config(config)
        else:
            raise ValueError(f"Unknown junction type: {config['junction_type']}")
        connections += [(f"V{vid}", junction.name) for vid in config["inlet_vessels"]]
        connections += [(junction.name, f"V{vid}") for vid in config["outlet_vessels"]]
        if junction.name in block_dict:
            raise RuntimeError(f"Junction {junction.name} already exists.")
        block_dict[junction.name] = junction

    return block_dict, connections


def create_vessel_blocks(parameters):
    """
    Purpose:
        Create the vessel blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
    Returns:
        void, but updates parameters["blocks"] to include the vessel_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    block_dict = {}
    connections = []
    for vessel_config in parameters["vessels"]:
        if vessel_config["zero_d_element_type"] == "BloodVessel":
            vessel = BloodVessel.from_config(vessel_config)
        else:
            raise NotImplementedError
        if "boundary_conditions" in vessel_config:
            if "inlet" in vessel_config["boundary_conditions"]:
                connections.append(
                    (
                        "BC" + str(vessel_config["vessel_id"]) + "_inlet",
                        vessel.name,
                    )
                )
            if "outlet" in vessel_config["boundary_conditions"]:
                connections.append(
                    (
                        vessel.name,
                        "BC" + str(vessel_config["vessel_id"]) + "_outlet",
                    )
                )
        if vessel.name in block_dict:
            raise RuntimeError(f"Vessel {vessel.name} already exists.")
        block_dict[vessel.name] = vessel
    return block_dict, connections


def create_bc_blocks(parameters, steady=False):
    """
    Purpose:
        Create the outlet bc (boundary condition) blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
    Returns:
        void, but updates parameters["blocks"] to include the blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    block_dict = {}
    for vessel_config in [
        v for v in parameters["vessels"] if "boundary_conditions" in v
    ]:
        locations = [
            loc
            for loc in ("inlet", "outlet")
            if loc in vessel_config["boundary_conditions"]
        ]
        for location in locations:
            vessel_id = vessel_config["vessel_id"]
            for config in parameters["boundary_conditions"]:
                if config["bc_name"] == vessel_config["boundary_conditions"][location]:
                    bc_config = dict(
                        name="BC" + str(vessel_id) + "_" + location,
                        steady=steady,
                        **config,
                    )
                    break

            if "t" in bc_config["bc_values"] and len(bc_config["bc_values"]["t"]) >= 2:
                cardiac_cycle_period = (
                    bc_config["bc_values"]["t"][-1] - bc_config["bc_values"]["t"][0]
                )
                if (
                    "cardiac_cycle_period" in parameters["simulation_parameters"]
                    and cardiac_cycle_period
                    != parameters["simulation_parameters"]["cardiac_cycle_period"]
                ):
                    raise RuntimeError(
                        f"The time series of the boundary condition for segment {vessel_id} does "
                        "not have the same cardiac cycle period as the other boundary conditions."
                    )
                elif "cardiac_cycle_period" not in parameters["simulation_parameters"]:
                    parameters["simulation_parameters"][
                        "cardiac_cycle_period"
                    ] = cardiac_cycle_period

            if bc_config["bc_type"] == "RESISTANCE":
                bc = ResistanceWithDistalPressure.from_config(bc_config)
            elif bc_config["bc_type"] == "RCR":
                bc = RCRBlockWithDistalPressure.from_config(bc_config)
            elif bc_config["bc_type"] == "FLOW":
                bc = FlowRef.from_config(bc_config)
            elif bc_config["bc_type"] == "PRESSURE":
                bc = PressureRef.from_config(bc_config)
            elif bc_config["bc_type"] == "CORONARY":
                bc = OpenLoopCoronaryWithDistalPressureBlock.from_config(bc_config)
            else:
                raise NotImplementedError
            block_dict[bc.name] = bc
    return block_dict


def create_blocks(parameters, steady=False):
    """
    Purpose:
        Create all LPNBlock objects for the 0d model.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
    Returns:
        void, but updates parameters to include:
            dict blocks
                = {block_name : block_object}
                -- where the values are the junction, vessel, outlet BC, and inlet BC block objects
    """
    junction_blocks, junction_connections = create_junction_blocks(parameters)
    vessel_blocks, vessel_connections = create_vessel_blocks(parameters)
    bc_blocks = create_bc_blocks(parameters, steady=steady)
    all_blocks = junction_blocks | vessel_blocks | bc_blocks
    dofhandler = DOFHandler()
    for ele1_name, ele2_name in junction_connections + vessel_connections:
        node = Node(
            all_blocks[ele1_name],
            all_blocks[ele2_name],
            name=ele1_name + "_" + ele2_name,
        )
        node.setup_dofs(dofhandler)
    for block in all_blocks.values():
        block.setup_dofs(dofhandler)
    return list(all_blocks.values()), dofhandler
