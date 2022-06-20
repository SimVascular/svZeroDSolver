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
from . import use_steady_bcs


class Wire:
    """
    Wires connect circuit elements and junctions
    They can only posses a single pressure and flow value (system variables)
    They can also only possess one element(or junction) at each end
    """

    def __init__(self, connecting_elements, name="NoNameWire"):
        self.name = name
        if len(connecting_elements) > 2:
            raise Exception(
                "Wire cannot connect to more than two elements at a time. Use a junction LPN block"
            )
        if not isinstance(connecting_elements, tuple):
            raise Exception("Connecting elements to wire should be passed as a 2-tuple")
        self.connecting_elements = connecting_elements
        self.LPN_solution_ids = [None] * 2


class LPNBlock:
    def __init__(self, connecting_block_list=None, name="NoName", flow_directions=[]):
        if connecting_block_list == None:
            connecting_block_list = []
        self.connecting_block_list = connecting_block_list
        self.num_connections = len(connecting_block_list)
        self.name = name
        self.neq = 2
        self.num_block_vars = 0

        self.inflow_wires = []
        self.outflow_wires = []

        # -1 : Inflow to block, +1 outflow from block
        self.flow_directions = flow_directions

        # solution IDs for the LPN block's internal solution variables
        self.LPN_solution_ids = []

        # block matrices
        self.mat = {}
        self.vec = {}

        # mat and vec assembly queue. To reduce the need to reassemble
        # matrices that havent't changed since last assembly, these attributes
        # are used to queue updated mats and vecs.
        self.mats_to_assemble = set()
        self.vecs_to_assemble = set()

        # row and column indices of block in global matrix
        self.global_col_id = []
        self.global_row_id = []

        self.flat_row_ids = []
        self.flat_col_ids = []

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

    def eqids(self, local_eq):
        # EqID returns variable's location in solution vector

        nwirevars = (
            self.num_connections * 2
        )  # num_connections is multipled by 2 because each wire has 2 soltns (P and Q)
        if local_eq < nwirevars:
            vtype = local_eq % 2  # 0 --> P, 1 --> Q
            wnum = int(local_eq / 2)

            # example: assume num_connections is 2. this can be a normal resistor block, which has 2 connections. then this R block has 2 connecting wires. thus, this R block has 4 related solution variables/unknowns (P_in, Q_in, P_out, Q_out). note that local_eq = local ID.
            #     then for these are the vtypes we get for each local_eq:
            #         local_eq    :     vtype     :     wnum
            #         0            :     0        :    0        <---    vtype = pressure, wnum = inlet wire
            #         1            :    1        :    0        <---    vtype = flow, wnum = inlet wire
            #         2            :    0        :    1        <---    vtype = pressure, wnum = outlet wire
            #         3            :    1        :    1        <---    vtype = flow, wnum = outlet wire
            #    note that vtype represents whether the solution variable in local_eq (local ID) is a P or Q solution
            #        and wnum represents whether the solution variable in local_eq comes from the inlet wire or the outlet wire, for this LPNBlock with 2 connections (one inlet, one outlet)
            return (self.outflow_wires + self.inflow_wires)[wnum].LPN_solution_ids[
                vtype
            ]
        else:  # this section will return the index at which the LPNBlock's  INTERNAL SOLUTION VARIABLES are stored in the global vector of solution unknowns/variables (i.e. I think RCR and OpenLoopCoronaryBlock have internal solution variables; these internal solution variables arent the P_in, Q_in, P_out, Q_out that correspond to the solutions on the attached wires, they are the solutions that are internal to the LPNBlock itself)
            vnum = local_eq - nwirevars
            return self.LPN_solution_ids[vnum]


class InternalJunction(LPNBlock):
    """
    Internal junction points between LPN blocks (for mesh refinement, does not appear as physical junction in model)
    """

    def __init__(
        self, connecting_block_list=None, name="NoNameJunction", flow_directions=None
    ):
        LPNBlock.__init__(
            self, connecting_block_list, name=name, flow_directions=flow_directions
        )
        self.neq = (
            self.num_connections
        )  # number of equations = num of blocks that connect to this junction, where the equations are 1) mass conservation 2) inlet pressures = outlet pressures

        # Number of variables per tuple = 2*num_connections
        # Number of equations = num_connections-1 Pressure equations, 1 flow equation
        # Format : P1,Q1,P2,Q2,P3,Q3, .., Pn,Qm
        self.mat["F"] = [
            (1.0,)
            + (0,) * (2 * i + 1)
            + (-1,)
            + (0,) * (2 * self.num_connections - 2 * i - 3)
            for i in range(self.num_connections - 1)
        ]

        tmp = (0,)
        for d in self.flow_directions[:-1]:
            tmp += (d,)
            tmp += (0,)

        tmp += (self.flow_directions[-1],)
        self.mat["F"].append(tmp)
        self.mat["F"] = np.array(self.mat["F"], dtype=float)
        self.mats_to_assemble.add("F")

    @classmethod
    def from_config(cls, config):
        name = config["junction_name"]
        if not name.startswith("J") and not name[1].isnumeric():
            raise ValueError(
                f"Invalid junction name {name}. Junction names must "
                "start with J following by a number."
            )
        connecting_block_list = []
        flow_directions = []
        for vessel_id in config["inlet_vessels"]:
            connecting_block_list.append("V" + str(vessel_id))
            flow_directions.append(-1)
        for vessel_id in config["outlet_vessels"]:
            connecting_block_list.append("V" + str(vessel_id))
            flow_directions.append(+1)
        if not (+1 in flow_directions) and (-1 in flow_directions):
            raise ValueError(
                "Junction block must have at least 1 inlet and 1 outlet " "connection."
            )
        return cls(connecting_block_list, name, flow_directions)


class BloodVesselJunction(InternalJunction):
    """
    Blood vessel junction (dummy for future implementation of blood pressure losses at junctions)
    """

    def __init__(
        self,
        j_params,
        connecting_block_list=None,
        name="NoNameJunction",
        flow_directions=None,
    ):
        InternalJunction.__init__(
            self, connecting_block_list, name=name, flow_directions=flow_directions
        )
        self.j_params = j_params


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
        connecting_block_list=None,
        name="NoNameBloodVessel",
        flow_directions=None,
    ):
        LPNBlock.__init__(
            self, connecting_block_list, name=name, flow_directions=flow_directions
        )
        self.neq = 3
        self.num_block_vars = 1
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
            connecting_block_list=config["connecting_blocks"],
            name=config["name"],
            flow_directions=config["flow_directions"],
        )

    def update_solution(self, sol):
        Q_in = np.abs(sol[self.inflow_wires[0].LPN_solution_ids[1]])
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
        connecting_block_list=None,
        name="NoNameUnsteadyResistanceWithDistalPressure",
        flow_directions=None,
    ):
        LPNBlock.__init__(
            self, connecting_block_list, name=name, flow_directions=flow_directions
        )
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
            connecting_block_list=config["connecting_blocks"],
            Rfunc=lambda t: config["bc_values"]["R"],
            Pref_func=lambda t: config["bc_values"]["Pd"],
            name=config["name"],
            flow_directions=config["flow_directions"],
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
        connecting_block_list=None,
        name="NoNameUnsteadyPressureRef",
        flow_directions=None,
    ):
        LPNBlock.__init__(
            self, connecting_block_list, name=name, flow_directions=flow_directions
        )
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
            connecting_block_list=config["connecting_blocks"],
            Pfunc=Pfunc,
            name=config["name"],
            flow_directions=config["flow_directions"],
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
        connecting_block_list=None,
        name="NoNameUnsteadyFlowRef",
        flow_directions=None,
    ):
        LPNBlock.__init__(
            self, connecting_block_list, name=name, flow_directions=flow_directions
        )
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
            connecting_block_list=config["connecting_blocks"],
            Qfunc=Qfunc,
            name=config["name"],
            flow_directions=config["flow_directions"],
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
        connecting_block_list=None,
        name="NoNameUnsteadyRCRBlockWithDistalPressure",
        flow_directions=None,
    ):
        LPNBlock.__init__(
            self, connecting_block_list, name=name, flow_directions=flow_directions
        )
        self.neq = 2
        self.num_block_vars = 1
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
            connecting_block_list=config["connecting_blocks"],
            name=config["name"],
            flow_directions=config["flow_directions"],
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
        connecting_block_list=None,
        name="NoNameCoronary",
        flow_directions=None,
    ):
        LPNBlock.__init__(
            self, connecting_block_list, name=name, flow_directions=flow_directions
        )
        self.neq = 2
        self.num_block_vars = 1
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
            connecting_block_list=config["connecting_blocks"],
            name=config["name"],
            flow_directions=config["flow_directions"],
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
    block_list = []
    for config in parameters["junctions"]:
        if config["junction_type"] in ["NORMAL_JUNCTION", "internal_junction"]:
            block_list.append(InternalJunction.from_config(config))
        elif config["junction_type"] == "BloodVesselJunction":
            block_list.append(BloodVesselJunction.from_config(config))
        else:
            raise ValueError(f"Unknown junction type: {config['junction_type']}")
    return block_list


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
    vessel_config = {}
    for vessel in parameters["vessels"]:
        vessel_id = vessel["vessel_id"]
        vessel_config[vessel_id] = dict(
            name="V" + str(vessel_id),
            connecting_blocks=[],
            flow_directions=[],
            **vessel,
        )
    for location in ["inlet", "outlet"]:
        ids_of_cap_vessels = use_steady_bcs.get_ids_of_cap_vessels(parameters, location)
        for vessel_id in ids_of_cap_vessels:
            flow_direction = -1 if location == "inlet" else +1
            vessel_config[vessel_id]["flow_directions"].append(flow_direction)
            vessel_config[vessel_id]["connecting_blocks"].append(
                f"BC{vessel_id}_{location}"
            )
        for junction in parameters["junctions"]:
            for vessel_id in junction[location + "_vessels"]:
                vessel_config[vessel_id]["connecting_blocks"].append(
                    junction["junction_name"]
                )
                flow_direction = +1 if location == "inlet" else -1
                vessel_config[vessel_id]["flow_directions"].append(flow_direction)
    block_list = []
    for vessel in vessel_config.values():
        if vessel["zero_d_element_type"] == "BloodVessel":
            block_list.append(BloodVessel.from_config(vessel))
        else:
            raise NotImplementedError
    return block_list


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
    block_list = []
    vessel_id_to_boundary_condition_map = {}
    for vessel in parameters["vessels"]:
        if "boundary_conditions" in vessel:
            vessel_id = vessel["vessel_id"]
            vessel_id_to_boundary_condition_map[vessel_id] = {}
            for location, bc_name in vessel["boundary_conditions"].items():
                for boundary_condition in parameters["boundary_conditions"]:
                    if boundary_condition["bc_name"] == bc_name:
                        vessel_id_to_boundary_condition_map[vessel_id][
                            location
                        ] = boundary_condition
    parameters[
        "vessel_id_to_boundary_condition_map"
    ] = vessel_id_to_boundary_condition_map
    for location in ["outlet", "inlet"]:
        for vessel_id in use_steady_bcs.get_ids_of_cap_vessels(parameters, location):
            block_name = "BC" + str(vessel_id) + "_" + location

            vessel_config = dict(
                name=block_name,
                connecting_blocks=["V" + str(vessel_id)],
                flow_directions=[-1] if location == "outlet" else [+1],
                **parameters["vessel_id_to_boundary_condition_map"][vessel_id][
                    location
                ],
            )

            if (
                "t" in vessel_config["bc_values"]
                and len(vessel_config["bc_values"]["t"]) >= 2
            ):
                cardiac_cycle_period = (
                    vessel_config["bc_values"]["t"][-1]
                    - vessel_config["bc_values"]["t"][0]
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

            if vessel_config["bc_type"] == "RESISTANCE":
                block_list.append(
                    UnsteadyResistanceWithDistalPressure.from_config(vessel_config)
                )
            elif vessel_config["bc_type"] == "RCR":
                block_list.append(
                    UnsteadyRCRBlockWithDistalPressure.from_config(vessel_config)
                )
            elif vessel_config["bc_type"] == "FLOW":
                block_list.append(UnsteadyFlowRef.from_config(vessel_config))
            elif vessel_config["bc_type"] == "PRESSURE":
                block_list.append(UnsteadyPressureRef.from_config(vessel_config))
            elif vessel_config["bc_type"] == "CORONARY":
                block_list.append(
                    OpenLoopCoronaryWithDistalPressureBlock.from_config(vessel_config)
                )

            else:
                raise NotImplementedError
    return block_list


def create_LPN_blocks(parameters):
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
    block_list = create_junction_blocks(parameters)
    block_list += create_vessel_blocks(parameters)
    block_list += create_bc_blocks(parameters)
    return block_list
