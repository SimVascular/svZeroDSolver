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


class wire:
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
        self.connecting_wires_list = []

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

    def add_connecting_wire(self, new_wire):
        self.connecting_wires_list.append(new_wire)

    def update_time(self, args):
        """
        Update time-dependent blocks
        """
        pass

    def update_solution(self, args):
        """
        Update solution-dependent blocks
        """
        pass

    def eqids(self, wire_dict, local_eq):
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

            return wire_dict[self.connecting_wires_list[wnum]].LPN_solution_ids[vtype]
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

    def update_solution(self, args):
        Q_in = np.abs(
            args["Solution"][
                args["Wire dictionary"][self.connecting_wires_list[0]].LPN_solution_ids[
                    1
                ]
            ]
        )
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

    def update_time(self, args):
        """
        the ordering is : (P_in,Q_in)
        """
        t = args["Time"]
        self.mat["F"][0, 1] = -self.Rfunc(t)
        self.vec["C"][0] = -self.Pref_func(t)
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

    def update_time(self, args):
        t = args["Time"]
        self.vec["C"][0] = -self.Pfunc(t)
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

    def update_time(self, args):
        t = args["Time"]
        self.vec["C"][0] = -self.Qfunc(t)
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

    def update_time(self, args):
        """
        unknowns = [P_in, Q_in, internal_var (Pressure at the intersection of the Rp, Rd, and C elements)]
        """
        t = args["Time"]
        Rd_t = self.Rd_func(t)
        self.mat["E"][1, 2] = -Rd_t * self.C_func(t)
        self.mat["F"][0, 1] = -self.Rp_func(t)
        self.mat["F"][1, 1] = Rd_t
        self.vec["C"][1] = self.Pref_func(t)
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

    def get_P_at_t(self, P, t):
        tt = P[:, 0]
        P_val = P[:, 1]
        _, td = divmod(t, self.cardiac_cycle_period)
        P_tt = np.interp(td, tt, P_val)
        return P_tt

    def update_time(self, args):
        # For this open-loop coronary BC, the ordering of solution unknowns is : (P_in, Q_in, V_im)
        # where V_im is the volume of the second capacitor, Cim
        # Q_in is the flow through the first resistor
        # and P_in is the pressure at the inlet of the first resistor
        ttt = args["Time"]
        Pim_value = self.get_P_at_t(self.Pim, ttt)
        Pv_value = self.get_P_at_t(self.Pv, ttt)
        self.vec["C"][0] = -self.Cim * Pim_value + self.Cim * Pv_value
        self.vec["C"][1] = (
            -self.Cim * (self.Rv + self.Ram) * Pim_value
            + self.Ram * self.Cim * Pv_value
        )
        self.vecs_to_assemble.add("C")
