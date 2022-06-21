# coding=utf-8

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

import numpy as np

from scipy.interpolate import CubicSpline


class DOFHandler:
    def __init__(self):
        self._var_counter = -1
        self._eqn_counter = -1
        self.variables = {}

    def register_variable(self, name=None):
        self._var_counter += 1
        if name is not None:
            self.variables[self._var_counter] = name
        return self._var_counter

    def register_equation(self):
        self._eqn_counter += 1
        return self._eqn_counter


class Wire:
    """
    Wires connect circuit elements and junctions
    They can only posses a single pressure and flow value (system variables)
    They can also only possess one element(or junction) at each end
    """

    def __init__(self, ele1, ele2, name="NoNameWire"):
        self.name = name
        self.flow_dof = None
        self.pres_dof = None
        ele1.outflow_wires.append(self)
        ele2.inflow_wires.append(self)

    def setup_dofs(self, dofhandler):
        self.flow_dof = dofhandler.register_variable("Q_" + self.name)
        self.pres_dof = dofhandler.register_variable("P_" + self.name)


class LPNBlock:
    def __init__(self, name="NoName"):

        self.name = name
        self.neq = 2

        self.inflow_wires = []
        self.outflow_wires = []

        # solution IDs for the LPN block's internal solution variables
        self.pres_dofs = []
        self.flow_dofs = []

        # block matrices
        self.mat = {}
        self.vec = {}

        # mat and vec assembly queue. To reduce the need to reassemble
        # matrices that havent't changed since last assembly, these attributes
        # are used to queue updated mats and vecs.
        self.mats_to_assemble = set()
        self.vecs_to_assemble = set()

        # row and column indices of block in global matrix
        self.global_col_id = None
        self.global_row_id = None

        self.flat_row_ids = None
        self.flat_col_ids = None

    def setup_dofs(self, dofhandler):
        self.flow_dofs = [
            dofhandler.register_variable("var_" + str(i) + "_" + self.name)
            for i in range(len(self.flow_dofs))
        ]

        self.pres_dofs = [
            dofhandler.register_variable("var_" + str(j) + "_" + self.name)
            for j in range(
                len(self.flow_dofs), len(self.pres_dofs) + len(self.flow_dofs)
            )
        ]

        self.global_col_id = []
        for wire in self.inflow_wires:
            self.global_col_id += [wire.pres_dof, wire.flow_dof]
        for wire in self.outflow_wires:
            self.global_col_id += [wire.pres_dof, wire.flow_dof]
        self.global_col_id += self.pres_dofs + self.flow_dofs

        self.global_row_id = [dofhandler.register_equation() for _ in range(self.neq)]

        meshgrid = np.array(
            np.meshgrid(self.global_row_id, self.global_col_id)
        ).T.reshape(-1, 2)
        self.flat_row_ids, self.flat_col_ids = meshgrid[:, 0], meshgrid[:, 1]

    def update_time(self, time):
        """
        Update time-dependent blocks
        """
        pass

    def update_solution(self, y):
        """
        Update solution-dependent blocks
        """
        pass


class InternalJunction(LPNBlock):
    """
    Internal junction points between LPN blocks (for mesh refinement, does not appear as physical junction in model)
    """

    def setup_dofs(self, dofhandler):
        num_inflow = len(self.inflow_wires)
        num_outflow = len(self.outflow_wires)
        self.neq = num_inflow + num_outflow
        super().setup_dofs(dofhandler)

        self.mat["F"] = np.zeros((self.neq, self.neq * 2))
        for i in range(self.neq - 1):
            self.mat["F"][i, [0, 2 * i + 2]] = [1.0, -1.0]
        self.mat["F"][-1, np.arange(1, 2 * self.neq, 2)] = [1.0] * num_inflow + [
            -1.0
        ] * num_outflow
        self.mats_to_assemble.add("F")

    @classmethod
    def from_config(cls, config):
        name = config["junction_name"]
        if not name.startswith("J") and not name[1].isnumeric():
            raise ValueError(
                f"Invalid junction name {name}. Junction names must "
                "start with J following by a number."
            )
        return cls(name)


class BloodVessel(LPNBlock):
    """
    Stenosis:
        equation: delta_P = ( K_t * rho / ( 2 * (A_0)**2 ) ) * ( ( A_0 / A_s ) - 1 )**2 * Q * abs(Q) + R_poiseuille * Q
                          =               stenosis_coefficient                          * Q * abs(Q) + R_poiseuille * Q

        source: Mirramezani, M., Shadden, S.C. A distributed lumped parameter model of blood flow. Annals of Biomedical Engineering. 2020.
    """

    def __init__(
        self,
        R,
        C,
        L,
        stenosis_coefficient,
        name="NoNameBloodVessel",
    ):
        LPNBlock.__init__(self, name=name)
        self.neq = 3
        self.pres_dofs = [None]
        self.R = R  # poiseuille resistance value = 8 * mu * L / (pi * r**4)
        self.C = C
        self.L = L
        self.stenosis_coefficient = stenosis_coefficient

        # the ordering of the solution variables is : (P_in, Q_in, P_out, Q_out)
        self.mat["E"] = np.zeros((3, 5), dtype=float)
        self.mat["E"][0, 3] = -self.L
        self.mat["E"][1, 4] = -self.C
        self.mat["F"] = np.array(
            [
                [1.0, 0.0, -1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, -1.0, 0.0],
                [1.0, 0.0, 0.0, 0.0, -1.0],
            ],
            dtype=float,
        )
        self.mat["dF"] = np.zeros((3, 5), dtype=float)

        # only necessary to assemble E in __init__, F and dF get assembled with update_solution
        self.mats_to_assemble.add("E")

    @classmethod
    def from_config(cls, config):
        return cls(
            R=config["zero_d_element_values"].get("R_poiseuille"),
            C=config["zero_d_element_values"].get("C", 0.0),
            L=config["zero_d_element_values"].get("L", 0.0),
            stenosis_coefficient=config["zero_d_element_values"].get(
                "stenosis_coefficient", 0.0
            ),
            name="V" + str(config["vessel_id"]),
        )

    def update_solution(self, sol):
        Q_in = np.abs(sol[self.inflow_wires[0].flow_dof])
        fac1 = -self.stenosis_coefficient * Q_in
        fac2 = fac1 - self.R
        self.mat["F"][[0, 2], 1] = fac2
        self.mat["dF"][[0, 2], 1] = fac1
        self.mats_to_assemble.update({"F", "dF"})


class UnsteadyResistanceWithDistalPressure(LPNBlock):
    def __init__(
        self,
        Rfunc,
        Pref_func,
        name="NoNameUnsteadyResistanceWithDistalPressure",
    ):
        LPNBlock.__init__(self, name=name)
        self.neq = 1
        self.Rfunc = Rfunc
        self.Pref_func = Pref_func
        self.mat["F"] = np.array(
            [
                [1.0, 0.0],
            ],
            dtype=float,
        )
        self.vec["C"] = np.array([0.0], dtype=float)

    @classmethod
    def from_config(cls, config):
        return cls(
            Rfunc=lambda t: config["bc_values"]["R"],
            Pref_func=lambda t: config["bc_values"]["Pd"],
            name=config["name"],
        )

    def update_time(self, time):
        """
        the ordering is : (P_in,Q_in)
        """
        self.mat["F"][0, 1] = -self.Rfunc(time)
        self.vec["C"][0] = -self.Pref_func(time)
        self.mats_to_assemble.add("F")
        self.vecs_to_assemble.add("C")


class UnsteadyPressureRef(LPNBlock):
    """
    Unsteady P reference
    """

    def __init__(
        self,
        Pfunc,
        name="NoNameUnsteadyPressureRef",
    ):
        LPNBlock.__init__(self, name=name)
        self.neq = 1
        self.Pfunc = Pfunc

        self.vec["C"] = np.zeros(1, dtype=float)
        self.mat["F"] = np.array([[1.0, 0.0]], dtype=float)

    @classmethod
    def from_config(cls, config):
        time = config["bc_values"]["t"]
        bc_values = config["bc_values"]["P"]
        Pfunc = CubicSpline(np.array(time), np.array(bc_values), bc_type="periodic")
        return cls(
            Pfunc=Pfunc,
            name=config["name"],
        )

    def update_time(self, time):
        self.vec["C"][0] = -self.Pfunc(time)
        self.vecs_to_assemble.add("C")
        self.mats_to_assemble.add("F")


class UnsteadyFlowRef(LPNBlock):
    """
    Flow reference
    """

    def __init__(
        self,
        Qfunc,
        name="NoNameUnsteadyFlowRef",
    ):
        LPNBlock.__init__(self, name=name)
        self.neq = 1
        self.Qfunc = Qfunc
        self.vec["C"] = np.zeros(1, dtype=float)
        self.mat["F"] = np.array([[0.0, 1.0]], dtype=float)

    @classmethod
    def from_config(cls, config):
        time = config["bc_values"]["t"]
        bc_values = config["bc_values"]["Q"]
        Qfunc = CubicSpline(np.array(time), np.array(bc_values), bc_type="periodic")
        return cls(
            Qfunc=Qfunc,
            name=config["name"],
        )

    def update_time(self, time):
        self.vec["C"][0] = -self.Qfunc(time)
        self.vecs_to_assemble.add("C")
        self.mats_to_assemble.add("F")


class UnsteadyRCRBlockWithDistalPressure(LPNBlock):
    """
    Unsteady RCR - time-varying RCR values
    Formulation includes additional variable : internal pressure proximal to capacitance.
    """

    def __init__(
        self,
        Rp_func,
        C_func,
        Rd_func,
        Pref_func,
        name="NoNameUnsteadyRCRBlockWithDistalPressure",
    ):
        LPNBlock.__init__(self, name=name)
        self.neq = 2
        self.pres_dofs = [None]
        self.Rp_func = Rp_func
        self.C_func = C_func
        self.Rd_func = Rd_func
        self.Pref_func = Pref_func

        self.mat["E"] = np.zeros((2, 3), dtype=float)
        self.mat["F"] = np.array([[1.0, 0.0, -1.0], [0.0, 0.0, -1.0]], dtype=float)
        self.vec["C"] = np.array([0.0, 0.0], dtype=float)

    @classmethod
    def from_config(cls, config):
        return cls(
            Rp_func=lambda t: config["bc_values"]["Rp"],
            C_func=lambda t: config["bc_values"]["C"],
            Rd_func=lambda t: config["bc_values"]["Rd"],
            Pref_func=lambda t: config["bc_values"]["Pd"],
            name=config["name"],
        )

    def update_time(self, time):
        """
        unknowns = [P_in, Q_in, internal_var (Pressure at the intersection of the Rp, Rd, and C elements)]
        """
        Rd_t = self.Rd_func(time)
        self.mat["E"][1, 2] = -Rd_t * self.C_func(time)
        self.mat["F"][0, 1] = -self.Rp_func(time)
        self.mat["F"][1, 1] = Rd_t
        self.vec["C"][1] = self.Pref_func(time)
        self.mats_to_assemble.update({"E", "F"})
        self.vecs_to_assemble.add("C")


class OpenLoopCoronaryWithDistalPressureBlock(LPNBlock):
    """
    open-loop coronary BC = RCRCR BC
    Publication reference: Kim, H. J. et al. Patient-specific modeling of blood flow and pressure in human coronary arteries. Annals of Biomedical Engineering 38, 3195–3209 (2010)."
    """

    def __init__(
        self,
        Ra,
        Ca,
        Ram,
        Cim,
        Rv,
        Pim,
        Pv,
        cardiac_cycle_period,
        name="NoNameCoronary",
    ):
        LPNBlock.__init__(self, name=name)
        self.neq = 2
        self.pres_dofs = [None]
        self.Ra = Ra
        self.Ca = Ca
        self.Ram = Ram
        self.Cim = Cim
        self.Rv = Rv
        self.Pim = Pim
        self.Pv = Pv
        self.cardiac_cycle_period = cardiac_cycle_period

        self.vec["C"] = np.zeros(2)
        self.mat["E"] = np.zeros((2, 3))
        self.mat["F"] = np.zeros((2, 3))
        self.mat["F"][0, 2] = -1.0

        Cim_Rv = self.Cim * self.Rv
        self.mat["E"][0, 0] = -self.Ca * Cim_Rv
        self.mat["E"][0, 1] = self.Ra * self.Ca * Cim_Rv
        self.mat["E"][0, 2] = -Cim_Rv
        self.mat["E"][1, 2] = -Cim_Rv * self.Ram
        self.mat["F"][0, 1] = Cim_Rv
        self.mat["F"][1, 0] = Cim_Rv
        self.mat["F"][1, 1] = -Cim_Rv * self.Ra
        self.mat["F"][1, 2] = -(self.Rv + self.Ram)
        self.mats_to_assemble.update({"E", "F"})

    @classmethod
    def from_config(cls, config):
        intramyocard_pres_time = config["bc_values"]["t"]
        bc_values_of_intramyocardial_pressure = config["bc_values"]["Pim"]
        Pim_func = np.zeros((len(intramyocard_pres_time), 2))
        Pv_distal_pressure_func = np.zeros((len(intramyocard_pres_time), 2))
        Pim_func[:, 0] = intramyocard_pres_time
        Pv_distal_pressure_func[:, 0] = intramyocard_pres_time
        Pim_func[:, 1] = bc_values_of_intramyocardial_pressure
        Pv_distal_pressure_func[:, 1] = (
            np.ones(len(intramyocard_pres_time)) * config["bc_values"]["P_v"]
        )

        return cls(
            Ra=config["bc_values"]["Ra1"],
            Ca=config["bc_values"]["Ca"],
            Ram=config["bc_values"]["Ra2"],
            Cim=config["bc_values"]["Cc"],
            Rv=config["bc_values"]["Rv1"],
            Pim=Pim_func,
            Pv=Pv_distal_pressure_func,
            cardiac_cycle_period=config["bc_values"]["t"][-1]
            - config["bc_values"]["t"][0],
            name=config["name"],
        )

    def get_P_at_t(self, P, t):
        tt = P[:, 0]
        P_val = P[:, 1]
        _, td = divmod(t, self.cardiac_cycle_period)
        P_tt = np.interp(td, tt, P_val)
        return P_tt

    def update_time(self, time):
        # For this open-loop coronary BC, the ordering of solution unknowns is : (P_in, Q_in, V_im)
        # where V_im is the volume of the second capacitor, Cim
        # Q_in is the flow through the first resistor
        # and P_in is the pressure at the inlet of the first resistor
        Pim_value = self.get_P_at_t(self.Pim, time)
        Pv_value = self.get_P_at_t(self.Pv, time)
        self.vec["C"][0] = -self.Cim * Pim_value + self.Cim * Pv_value
        self.vec["C"][1] = (
            -self.Cim * (self.Rv + self.Ram) * Pim_value
            + self.Ram * self.Cim * Pv_value
        )
        self.vecs_to_assemble.add("C")


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
            junction = InternalJunction.from_config(config)
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


def create_bc_blocks(parameters):
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
                        name="BC" + str(vessel_id) + "_" + location, **config
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
                bc = UnsteadyResistanceWithDistalPressure.from_config(bc_config)
            elif bc_config["bc_type"] == "RCR":
                bc = UnsteadyRCRBlockWithDistalPressure.from_config(bc_config)
            elif bc_config["bc_type"] == "FLOW":
                bc = UnsteadyFlowRef.from_config(bc_config)
            elif bc_config["bc_type"] == "PRESSURE":
                bc = UnsteadyPressureRef.from_config(bc_config)
            elif bc_config["bc_type"] == "CORONARY":
                bc = OpenLoopCoronaryWithDistalPressureBlock.from_config(bc_config)
            else:
                raise NotImplementedError
            block_dict[bc.name] = bc
    return block_dict


def create_blocks(parameters):
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
    bc_blocks = create_bc_blocks(parameters)
    all_blocks = junction_blocks | vessel_blocks | bc_blocks
    dofhandler = DOFHandler()
    for ele1_name, ele2_name in junction_connections + vessel_connections:
        wire = Wire(
            all_blocks[ele1_name],
            all_blocks[ele2_name],
            name=ele1_name + "_" + ele2_name,
        )
        wire.setup_dofs(dofhandler)
    for block in all_blocks.values():
        block.setup_dofs(dofhandler)
    return list(all_blocks.values()), dofhandler
