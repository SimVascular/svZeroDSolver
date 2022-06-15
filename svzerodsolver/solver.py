#!/usr/bin/env python3
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

"""
This code simulates the 0D model.

Available vessel modeling types:
    1) R (constant resistor)
    2) C (constant capacitor)
    3) L (constant inductor)
    4) RC (constant resistor-capacitor)
    5) RL (constant resistor-inductor)
    6) RCL (constant resistor-capacitor-inductor)
    7) custom, user-defined vessel model

Available junction modeling types:
    1) NORMAL_JUNCTION (mass conservation and pressure continuity only)
    2) custom, user-defined junction type

Available boundary conditions (BCs):
- inlet BCs:
        1) constant flow rate
        2) time-dependent flow rate
        3) constant pressure
        4) time-dependent pressure
        5) custom, user-defined inlet BC
- outlet BCs:
        1) constant resistor with user-prescribed constant distal pressure value
        2) time-dependent resistor with user-prescribed time-dependent distal pressure value
        3) constant RCR with a distal pressure of zero
        4) time-dependent RCR with a distal pressure of zero
        5) constant flow rate
        6) time-dependent flow rate
        7) constant pressure
        8) time-dependent pressure
        9) steady coronary with time-dependent intramyocadial pressure with user-prescribed constant distal pressure
        10) custom, user-defined outlet BC
"""
import sys
import re
import os
import copy
import numpy as np
import scipy.interpolate
import json

from . import blocks as ntwku
from . import connections
from . import time_integration as time_int
from . import use_steady_bcs

import argparse
from collections import defaultdict


def create_unsteady_bc_value_function(time, bc_values):
    """
    Purpose:
        Create a time-dependent function for bc_values (boundary condition values) by interpolating bc_values over time with a cubic spline and periodic boundaries
    Caveats:
        1) bc_values must have the same start and end values, since periodic boundariess are used for the cubic spline
        2) bc_values and time must have the same length
    Inputs:
        list time
            = [list of time points]
        list bc_values
            = [list of bc_values]
    Returns:
        cubic spline for bc_values as a function of time
    Reference:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html
    """
    if len(bc_values) == 2 and bc_values[0] == bc_values[1]:

        def function(t):
            """
            Constant boundary condition. Need to provide a "t" (time) input to work with the general (unsteady) blocks defined in blocks.py
            """
            return bc_values[0]

        return function
    else:
        return scipy.interpolate.CubicSpline(
            np.array(time), np.array(bc_values), bc_type="periodic"
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
    junction_blocks = {}  # {block_name : block_object}
    for junction in parameters["junctions"]:
        junction_name = junction["junction_name"]
        if not junction_name.startswith("J") and not junction_name[1].isnumeric():
            message = (
                "Error. Joint name, "
                + junction_name
                + ", is not 'J' followed by numeric values. The 0D solver assumes that all joint names are 'J' followed by numeric values in the 0d solver input file. Note that the joint names are the same as the junction names."
            )
            raise RuntimeError(message)
        connecting_block_list = []
        flow_directions = []
        for vessel_id in junction["inlet_vessels"]:
            connecting_block_list.append("V" + str(vessel_id))
            flow_directions.append(-1)
        for vessel_id in junction["outlet_vessels"]:
            connecting_block_list.append("V" + str(vessel_id))
            flow_directions.append(+1)
        if (+1 in flow_directions) and (-1 in flow_directions):
            if junction["junction_type"] in ["NORMAL_JUNCTION", "internal_junction"]:
                junction_blocks[junction_name] = ntwku.InternalJunction(
                    connecting_block_list=connecting_block_list,
                    name=junction_name,
                    flow_directions=flow_directions,
                )
            elif junction["junction_type"] == "BloodVesselJunction":
                junction_blocks[junction_name] = ntwku.BloodVesselJunction(
                    junction,
                    connecting_block_list=connecting_block_list,
                    name=junction_name,
                    flow_directions=flow_directions,
                )
            else:
                raise ValueError("Unknown junction type " + junction["junction_type"])
        else:
            message = (
                "Error. Junction block, "
                + junction_name
                + ", must have at least 1 inlet connection and 1 outlet connection."
            )
            raise RuntimeError(message)
    parameters["blocks"].update(junction_blocks)


def get_vessel_ids_list(parameters):
    vessel_ids_list = []
    for vessel in parameters["vessels"]:
        vessel_ids_list.append(vessel["vessel_id"])
    return vessel_ids_list


def get_vessel_block_helpers(parameters):
    """
    Purpose:
        Create helper dictionaries to support the creation of the vessel blocks.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
    Returns:
        dict vessel_blocks_connecting_block_lists
            = {vessel_id : connecting_block_list}
        dict vessel_blocks_flow_directions
            = {vessel_id : flow_directions}
        dict vessel_blocks_names
            = {vessel_id : block_name}
    """
    vessel_blocks_connecting_block_lists = {}
    vessel_blocks_flow_directions = {}
    vessel_blocks_names = {}
    vessel_ids_list = get_vessel_ids_list(parameters)
    for vessel_id in vessel_ids_list:
        vessel_blocks_connecting_block_lists[vessel_id] = []
        vessel_blocks_flow_directions[vessel_id] = []
        vessel_blocks_names[vessel_id] = "V" + str(vessel_id)
    for location in ["inlet", "outlet"]:
        ids_of_cap_vessels = use_steady_bcs.get_ids_of_cap_vessels(parameters, location)
        for vessel_id in ids_of_cap_vessels:
            if location == "inlet":
                bc_block_name = "BC" + str(vessel_id) + "_inlet"
                vessel_blocks_flow_directions[vessel_id].append(-1)
            else:
                bc_block_name = "BC" + str(vessel_id) + "_outlet"
                vessel_blocks_flow_directions[vessel_id].append(+1)
            vessel_blocks_connecting_block_lists[vessel_id].append(bc_block_name)
    for junction in parameters["junctions"]:
        junction_name = junction["junction_name"]
        if not junction_name.startswith("J") and not junction_name[1].isnumeric():
            message = (
                "Joint name, "
                + junction_name
                + ", is not 'J' followed by numeric values. The 0D solver assumes that all joint names are 'J' followed by numeric values in the 0d solver input file. Note that the joint names are the same as the junction names."
            )
            raise RuntimeError(message)
        for location in ["inlet", "outlet"]:
            for vessel_id in junction[location + "_vessels"]:
                vessel_blocks_connecting_block_lists[vessel_id].append(junction_name)
                if location == "inlet":
                    vessel_blocks_flow_directions[vessel_id].append(+1)
                else:
                    vessel_blocks_flow_directions[vessel_id].append(-1)
    return (
        vessel_blocks_connecting_block_lists,
        vessel_blocks_flow_directions,
        vessel_blocks_names,
    )


def create_vessel_id_to_zero_d_element_dict(parameters):
    vessel_id_to_zero_d_element_dict = {}
    for vessel in parameters["vessels"]:
        vessel_id = vessel["vessel_id"]
        zero_d_element_type = vessel["zero_d_element_type"]
        zero_d_element_values = vessel["zero_d_element_values"]
        vessel_id_to_zero_d_element_dict[vessel_id] = {
            "zero_d_element_type": zero_d_element_type,
            "zero_d_element_values": zero_d_element_values,
        }
    return vessel_id_to_zero_d_element_dict


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
    vessel_blocks = {}  # {block_name : block_object}
    (
        vessel_blocks_connecting_block_lists,
        vessel_blocks_flow_directions,
        vessel_blocks_names,
    ) = get_vessel_block_helpers(parameters)
    vessel_ids_list = get_vessel_ids_list(parameters)
    vessel_id_to_zero_d_element_dict = create_vessel_id_to_zero_d_element_dict(
        parameters
    )
    for vessel_id in vessel_ids_list:
        block_name = vessel_blocks_names[vessel_id]
        connecting_block_list = vessel_blocks_connecting_block_lists[vessel_id]
        flow_directions = vessel_blocks_flow_directions[vessel_id]
        zero_d_element_type = vessel_id_to_zero_d_element_dict[vessel_id][
            "zero_d_element_type"
        ]
        if zero_d_element_type == "BloodVessel":
            if (
                "R_poiseuille"
                not in vessel_id_to_zero_d_element_dict[vessel_id][
                    "zero_d_element_values"
                ]
            ):
                message = "Error. BloodVessel block requires viscous Poiseuille-based resistance."
                raise RuntimeError(message)
            R = vessel_id_to_zero_d_element_dict[vessel_id]["zero_d_element_values"][
                "R_poiseuille"
            ]
            C = (
                vessel_id_to_zero_d_element_dict[vessel_id]["zero_d_element_values"][
                    "C"
                ]
                if "C"
                in vessel_id_to_zero_d_element_dict[vessel_id]["zero_d_element_values"]
                else 0
            )
            L = (
                vessel_id_to_zero_d_element_dict[vessel_id]["zero_d_element_values"][
                    "L"
                ]
                if "L"
                in vessel_id_to_zero_d_element_dict[vessel_id]["zero_d_element_values"]
                else 0
            )
            stenosis_coefficient = (
                vessel_id_to_zero_d_element_dict[vessel_id]["zero_d_element_values"][
                    "stenosis_coefficient"
                ]
                if "stenosis_coefficient"
                in vessel_id_to_zero_d_element_dict[vessel_id]["zero_d_element_values"]
                else 0
            )
            vessel_blocks[block_name] = ntwku.BloodVessel(
                R=R,
                C=C,
                L=L,
                stenosis_coefficient=stenosis_coefficient,
                connecting_block_list=connecting_block_list,
                name=block_name,
                flow_directions=flow_directions,
            )
        else:  # this is a custom, user-defined element block
            raise NotImplementedError
    parameters["blocks"].update(vessel_blocks)


def create_outlet_bc_blocks(parameters):
    """
    Purpose:
        Create the outlet bc (boundary condition) blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
    Returns:
        void, but updates parameters["blocks"] to include the outlet_bc_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    outlet_bc_blocks = {}  # {block_name : block_object}
    outlet_vessels_of_model = use_steady_bcs.get_ids_of_cap_vessels(
        parameters, "outlet"
    )
    vessel_id_to_boundary_condition_map = parameters[
        "vessel_id_to_boundary_condition_map"
    ]
    for vessel_id in outlet_vessels_of_model:
        block_name = "BC" + str(vessel_id) + "_outlet"
        connecting_block_list = ["V" + str(vessel_id)]
        flow_directions = [-1]

        if (
            vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_type"]
            == "RESISTANCE"
        ):
            R = vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_values"][
                "R"
            ]
            R_func = create_unsteady_bc_value_function([0.0, 1.0], [R, R])

            Pref = vessel_id_to_boundary_condition_map[vessel_id]["outlet"][
                "bc_values"
            ]["Pd"]

            Pref_func = create_unsteady_bc_value_function([0.0, 1.0], [Pref, Pref])
            outlet_bc_blocks[block_name] = ntwku.UnsteadyResistanceWithDistalPressure(
                connecting_block_list=connecting_block_list,
                Rfunc=R_func,
                Pref_func=Pref_func,
                name=block_name,
                flow_directions=flow_directions,
            )
        elif (
            vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_type"] == "RCR"
        ):
            Rp = vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_values"][
                "Rp"
            ]
            Rp_func = create_unsteady_bc_value_function([0.0, 1.0], [Rp, Rp])

            C = vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_values"][
                "C"
            ]
            C_func = create_unsteady_bc_value_function([0.0, 1.0], [C, C])

            Rd = vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_values"][
                "Rd"
            ]
            Rd_func = create_unsteady_bc_value_function([0.0, 1.0], [Rd, Rd])

            Pref = vessel_id_to_boundary_condition_map[vessel_id]["outlet"][
                "bc_values"
            ]["Pd"]
            Pref_func = create_unsteady_bc_value_function([0.0, 1.0], [Pref, Pref])

            outlet_bc_blocks[block_name] = ntwku.UnsteadyRCRBlockWithDistalPressure(
                Rp_func=Rp_func,
                C_func=C_func,
                Rd_func=Rd_func,
                Pref_func=Pref_func,
                connecting_block_list=connecting_block_list,
                name=block_name,
                flow_directions=flow_directions,
            )
        elif (
            vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_type"]
            == "FLOW"
        ):
            time = vessel_id_to_boundary_condition_map[vessel_id]["outlet"][
                "bc_values"
            ]["t"]
            bc_values = vessel_id_to_boundary_condition_map[vessel_id]["outlet"][
                "bc_values"
            ]["Q"]
            Qfunc = create_unsteady_bc_value_function(time, bc_values)
            outlet_bc_blocks[block_name] = ntwku.UnsteadyFlowRef(
                connecting_block_list=connecting_block_list,
                Qfunc=Qfunc,
                name=block_name,
                flow_directions=flow_directions,
            )
        elif (
            vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_type"]
            == "PRESSURE"
        ):
            time = vessel_id_to_boundary_condition_map[vessel_id]["outlet"][
                "bc_values"
            ]["t"]
            bc_values = vessel_id_to_boundary_condition_map[vessel_id]["outlet"][
                "bc_values"
            ]["P"]
            Pfunc = create_unsteady_bc_value_function(time, bc_values)
            outlet_bc_blocks[block_name] = ntwku.UnsteadyPressureRef(
                connecting_block_list=connecting_block_list,
                Pfunc=Pfunc,
                name=block_name,
                flow_directions=flow_directions,
            )
        elif (
            vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_type"]
            == "CORONARY"
        ):
            "Publication reference: Kim, H. J. et al. Patient-specific modeling of blood flow and pressure in human coronary arteries. Annals of Biomedical Engineering 38, 3195–3209 (2010)."
            Ra1 = vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_values"][
                "Ra1"
            ]
            Ra2 = vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_values"][
                "Ra2"
            ]
            Ca = vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_values"][
                "Ca"
            ]
            Cc = vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_values"][
                "Cc"
            ]
            Rv1 = vessel_id_to_boundary_condition_map[vessel_id]["outlet"]["bc_values"][
                "Rv1"
            ]
            Pv_distal_pressure = vessel_id_to_boundary_condition_map[vessel_id][
                "outlet"
            ]["bc_values"]["P_v"]

            time_of_intramyocardial_pressure = vessel_id_to_boundary_condition_map[
                vessel_id
            ]["outlet"]["bc_values"]["t"]

            bc_values_of_intramyocardial_pressure = vessel_id_to_boundary_condition_map[
                vessel_id
            ]["outlet"]["bc_values"]["Pim"]

            if "cardiac_cycle_period" in parameters["simulation_parameters"]:
                if (
                    time_of_intramyocardial_pressure[-1]
                    - time_of_intramyocardial_pressure[0]
                    != parameters["simulation_parameters"]["cardiac_cycle_period"]
                ):
                    message = (
                        "Error. The time series of the intramyocadial pressure for the coronary boundary condition for segment #"
                        + str(vessel_id)
                        + " does not have the same cardiac cycle period as the other boundary conditions.  All boundary conditions, including the inlet and outlet boundary conditions, should have the same prescribed cardiac cycle period. Note that each boundary conditions must be prescribed over exactly one cardiac cycle."
                    )  # todo: fix bug where this code does not work if the user prescribes only a single time point for the intramyocadial pressure time history (the code should work though b/c prescription of a single time point for the intramyocardial pressure suggests a steady intramyocadial pressure)
                    raise RuntimeError(message)
            else:
                parameters["simulation_parameters"].update(
                    {
                        "cardiac_cycle_period": time_of_intramyocardial_pressure[-1]
                        - time_of_intramyocardial_pressure[0]
                    }
                )

            Pim_func = np.zeros((len(time_of_intramyocardial_pressure), 2))
            Pv_distal_pressure_func = np.zeros(
                (len(time_of_intramyocardial_pressure), 2)
            )

            Pim_func[:, 0] = time_of_intramyocardial_pressure
            Pv_distal_pressure_func[:, 0] = time_of_intramyocardial_pressure

            Pim_func[:, 1] = bc_values_of_intramyocardial_pressure
            Pv_distal_pressure_func[:, 1] = (
                np.ones(len(time_of_intramyocardial_pressure)) * Pv_distal_pressure
            )

            outlet_bc_blocks[
                block_name
            ] = ntwku.OpenLoopCoronaryWithDistalPressureBlock(
                Ra=Ra1,
                Ca=Ca,
                Ram=Ra2,
                Cim=Cc,
                Rv=Rv1,
                Pim=Pim_func,
                Pv=Pv_distal_pressure_func,
                cardiac_cycle_period=parameters["simulation_parameters"][
                    "cardiac_cycle_period"
                ],
                connecting_block_list=connecting_block_list,
                name=block_name,
                flow_directions=flow_directions,
            )

        else:  # this is a custom, user-defined outlet bc block
            raise NotImplementedError
    parameters["blocks"].update(outlet_bc_blocks)


def create_inlet_bc_blocks(parameters):
    """
    Purpose:
        Create the inlet bc (boundary condition) blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
    Returns:
        void, but updates parameters["blocks"] to include the inlet_bc_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    inlet_bc_blocks = {}  # {block_name : block_object}
    inlet_vessels_of_model = use_steady_bcs.get_ids_of_cap_vessels(parameters, "inlet")
    vessel_id_to_boundary_condition_map = parameters[
        "vessel_id_to_boundary_condition_map"
    ]
    for vessel_id in inlet_vessels_of_model:
        block_name = "BC" + str(vessel_id) + "_inlet"
        connecting_block_list = ["V" + str(vessel_id)]
        flow_directions = [+1]
        if vessel_id_to_boundary_condition_map[vessel_id]["inlet"]["bc_type"] == "FLOW":
            time = vessel_id_to_boundary_condition_map[vessel_id]["inlet"]["bc_values"][
                "t"
            ]
            bc_values = vessel_id_to_boundary_condition_map[vessel_id]["inlet"][
                "bc_values"
            ]["Q"]
            Qfunc = create_unsteady_bc_value_function(time, bc_values)
            inlet_bc_blocks[block_name] = ntwku.UnsteadyFlowRef(
                Qfunc=Qfunc,
                connecting_block_list=connecting_block_list,
                name=block_name,
                flow_directions=flow_directions,
            )
            if len(time) >= 2:
                cardiac_cycle_period = time[-1] - time[0]
            if "cardiac_cycle_period" in parameters["simulation_parameters"]:
                if (
                    cardiac_cycle_period
                    != parameters["simulation_parameters"]["cardiac_cycle_period"]
                ):
                    message = (
                        "Error. The time series of the boundary condition for segment #"
                        + str(vessel_id)
                        + " does not have the same cardiac cycle period as the other boundary conditions.  All boundary conditions, including the inlet and outlet boundary conditions, should have the same prescribed cardiac cycle period. Note that each boundary conditions must be prescribed over exactly one cardiac cycle."
                    )
                    raise RuntimeError(message)
            else:
                parameters["simulation_parameters"].update(
                    {"cardiac_cycle_period": cardiac_cycle_period}
                )
        elif (
            vessel_id_to_boundary_condition_map[vessel_id]["inlet"]["bc_type"]
            == "PRESSURE"
        ):
            time = vessel_id_to_boundary_condition_map[vessel_id]["inlet"]["bc_values"][
                "t"
            ]
            bc_values = vessel_id_to_boundary_condition_map[vessel_id]["inlet"][
                "bc_values"
            ]["P"]
            Pfunc = create_unsteady_bc_value_function(time, bc_values)
            inlet_bc_blocks[block_name] = ntwku.UnsteadyPressureRef(
                Pfunc=Pfunc,
                connecting_block_list=connecting_block_list,
                name=block_name,
                flow_directions=flow_directions,
            )
            if len(time) >= 2:
                cardiac_cycle_period = time[-1] - time[0]
            if "cardiac_cycle_period" in parameters["simulation_parameters"]:
                if (
                    cardiac_cycle_period
                    != parameters["simulation_parameters"]["cardiac_cycle_period"]
                ):
                    message = (
                        "Error. The time series of the boundary condition for segment #"
                        + str(vessel_id)
                        + " does not have the same cardiac cycle period as the other boundary conditions.  All boundary conditions, including the inlet and outlet boundary conditions, should have the same prescribed cardiac cycle period. Note that each boundary conditions must be prescribed over exactly one cardiac cycle."
                    )
                    raise RuntimeError(message)
            else:
                parameters["simulation_parameters"].update(
                    {"cardiac_cycle_period": cardiac_cycle_period}
                )
        else:  # this is a custom, user-defined inlet bc block
            raise NotImplementedError
    parameters["blocks"].update(inlet_bc_blocks)


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
    blocks = {}  # {block_name : block_object}
    parameters.update({"blocks": blocks})
    create_junction_blocks(parameters)
    create_vessel_blocks(parameters)
    use_steady_bcs.create_vessel_id_to_boundary_condition_map(parameters)
    create_outlet_bc_blocks(parameters)
    create_inlet_bc_blocks(parameters)
    parameters.update({"block_names": list(parameters["blocks"].keys())})


def set_solver_parameters(parameters):
    """
    Purpose:
        Set the 0d simulation time-stepping parameters
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
    Returns:
        void, but updates parameters to include:
            float delta_t
                = constant time step size for the 0d simulation
            int total_number_of_simulated_time_steps
                = total number of time steps to simulate for the entire 0d simulation
    """
    if "cardiac_cycle_period" not in parameters["simulation_parameters"]:
        parameters["simulation_parameters"].update(
            {"cardiac_cycle_period": 1.0}
        )  # default period of cardiac cycle [sec]
    delta_t = parameters["simulation_parameters"]["cardiac_cycle_period"] / (
        parameters["simulation_parameters"]["number_of_time_pts_per_cardiac_cycle"] - 1
    )
    parameters["simulation_parameters"].update({"delta_t": delta_t})
    total_number_of_simulated_time_steps = int(
        (
            parameters["simulation_parameters"]["number_of_time_pts_per_cardiac_cycle"]
            - 1
        )
        * parameters["simulation_parameters"]["number_of_cardiac_cycles"]
        + 1
    )
    parameters["simulation_parameters"].update(
        {"total_number_of_simulated_time_steps": total_number_of_simulated_time_steps}
    )


def load_in_ics(var_name_list, ICs_dict):

    var_name_list_loaded = ICs_dict["var_name_list"]
    y_initial_loaded = ICs_dict["y"]
    ydot_initial_loaded = ICs_dict["ydot"]

    y_initial = np.zeros(len(var_name_list))
    ydot_initial = np.zeros(len(var_name_list))
    for i in range(len(var_name_list)):
        var_name = var_name_list[i]
        ind = var_name_list_loaded.index(var_name)
        y_initial[i] = y_initial_loaded[ind]
        ydot_initial[i] = ydot_initial_loaded[ind]

    return y_initial, ydot_initial


def run_network_util(
    parameters,
    use_ICs_from_npy_file,
    ICs_npy_file_path,
    save_y_ydot_to_npy,
    y_ydot_file_path,
    simulation_start_time,
):
    """
    Purpose:
        Run functions from network_util_NR to execute the 0d simulation and generate simulation results (pressure, flow rates).
    Inputs:
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
        boolean draw_directed_graph
            = True to visualize the 0d model as a directed graph using networkx -- saves the graph to a .png file (hierarchical graph layout) and a networkx .dot file; False, otherwise. .dot file can be opened with neato from graphviz to visualize the directed in a different format.
        boolean use_ICs_from_npy_file
            = True to use user-prescribed ICs (saved in a .npy file)
        string ICs_npy_file_path
            = path to .npy file storing the initial conditions
        boolean save_y_ydot_to_npy
            = True to save the entire 0d solution and its time derivative at the final time step to a .npy file
        string y_ydot_file_path
            = name of the .npy file to which the 0d solution and its time derivative at the final time step will be saved
        float simulation_start_time
            = initial time of 0d simulation
    Returns:
        np.array zero_d_time
            = np.array of simulated time points
        np.array results_0d
            = np.array of the 0d simulation results, where the rows correspond to the each time point in zero_d_time and the column j corresponds to the solution for item j in var_name_list
        list var_name_list
            = list of the names of the 0d simulation results; most of the items in var_name_list are the QoIs + the names of the wires used in the 0d model (the wires connect the 0d blocks), where the wire names are usually comprised of the wire's inlet block name + "_" + the wire's outlet block name
                Example:
                    for var_name_list = ['P_V6_BC6_outlet', 'Q_V6_BC6_outlet'], then results_0d[:, i] holds the pressure (i = 0) or flow rate simulation result (i = 1) (both as np.arrays) for wire R6_BC6, which corresponds to cap segment #6
        np.array y_next
            = 0d solution at the final time step
        np.array ydot_next
            = time derivative of the 0d solution at the final time step
        var_name_list_original
            = list of the names of the 0d simulation results; most of the items in var_name_list are the QoIs + the names of the wires used in the 0d model (the wires connect the 0d blocks), where the wire names are usually comprised of the wire's inlet block name + "_" + the wire's outlet block name
    """

    block_list = list(parameters["blocks"].values())
    wire_dict = connections.connect_blocks_by_inblock_list(block_list)
    neq = connections.compute_neq(
        block_list, wire_dict
    )  # number of equations governing the 0d model
    for block in block_list:  # run a consistency check
        connections.check_block_connection(block)
    var_name_list = connections.assign_global_ids(
        block_list, wire_dict
    )  # assign solution variables with global ID

    # initialize solution structures
    if use_ICs_from_npy_file:
        ICs_dict = np.load(ICs_npy_file_path, allow_pickle=True).item()
        y_initial, ydot_initial = load_in_ics(var_name_list, ICs_dict)
    else:
        y_initial, ydot_initial = connections.initialize_solution_structures(
            neq
        )  # initial conditions for all solutions are zero
    y_next = y_initial.copy()
    ydot_next = ydot_initial.copy()

    rho = 0.1

    print("starting simulation")

    ylist = [y_next.copy()]
    parameters["initial_time"] = simulation_start_time
    dt = parameters["simulation_parameters"]["delta_t"]
    tlist = np.arange(
        parameters["initial_time"],
        parameters["simulation_parameters"]["total_number_of_simulated_time_steps"]
        * dt,
        dt,
    )

    # create time integration
    t_int = time_int.GeneralizedAlpha(
        rho, y_next.shape[0], parameters["simulation_parameters"]["delta_t"]
    )

    for time in tlist[:-1]:
        y_next, ydot_next = t_int.step(
            y_next,
            ydot_next,
            time,
            block_list,
        )
        ylist.append(y_next)

    if save_y_ydot_to_npy:
        save_ics(y_ydot_file_path, y_next, ydot_next, var_name_list)

    var_name_list_original = copy.deepcopy(var_name_list)
    results_0d = np.array(ylist)
    return (
        tlist,
        results_0d,
        var_name_list,
        y_next,
        ydot_next,
        var_name_list_original,
    )


def save_ics(y_ydot_file_path, y_next, ydot_next, var_name_list):
    np.save(
        y_ydot_file_path,
        {"y": y_next, "ydot": ydot_next, "var_name_list": var_name_list},
    )


def reformat_network_util_results_all(zero_d_time, results_0d, var_name_list):
    """
    Purpose:
        Reformat all 0d simulation results (results_0d) into a dictionary (zero_d_results_for_var_names)
    Inputs:
        np.array zero_d_time
            = np.array of simulated time points
        np.array results_0d
            = np.array of the 0d simulation results, where the rows correspond to the each simulated time point and column j corresponds to the 0d solution for the solution variable name in var_name_list[j]
        list var_name_list
            = list of the 0d simulation results' solution variable names; most of the items in var_name_list are the QoIs + the names of the wires used in the 0d model (the wires connecting the 0d LPNBlock objects), where the wire names are usually comprised of the wire's inlet block name + "_" + the wire's outlet block name
                Example:
                    for var_name_list = ['P_V6_BC6_outlet', 'Q_V6_BC6_outlet'], then results_0d[:, i] holds the pressure (i = 0) or flow rate simulation result (i = 1) (both as np.arrays) for wire R6_BC6_outlet. This wire connects a resistance vessel block to an outlet BC block (specifically for vessel segment #6)
    Returns:
        dict zero_d_results_for_var_names
            =   {
                    "time" : np.array of simulated time points,

                    "flow" : {var_name : np.array of flow rate,

                    "pressure" : {var_name : np.array of pressure},

                    "internal" : {var_name : np.array of internal block solutions},

                        where var_name is an item in var_name_list (var_name_list generated from run_network_util)
                }
    """
    zero_d_results_for_var_names = {
        "flow": {},
        "pressure": {},
        "time": zero_d_time,
        "internal": {},
    }
    for i in range(len(var_name_list)):
        var_name = var_name_list[i]
        res = results_0d[:, i]
        if var_name.startswith("Q_"):
            zero_d_results_for_var_names["flow"][var_name] = res
        elif var_name.startswith("P_"):
            zero_d_results_for_var_names["pressure"][var_name] = res
        elif var_name.startswith("var_"):
            zero_d_results_for_var_names["internal"][var_name] = res
        else:
            message = (
                "Error. There are unaccounted for solution variables here, for var_name = "
                + var_name
            )
            raise RuntimeError(message)
    return zero_d_results_for_var_names


def get_vessel_id_to_length_map(parameters):
    vessel_id_to_length_map = {}
    for vessel in parameters["vessels"]:
        vessel_id = vessel["vessel_id"]
        vessel_length = vessel["vessel_length"]
        vessel_id_to_length_map[vessel_id] = vessel_length
    return vessel_id_to_length_map


def initialize_0d_results_dict_branch(parameters, zero_d_time):
    zero_d_results = {"flow": {}, "pressure": {}, "distance": {}}
    num_time_pts = len(zero_d_time)
    branch_segment_ids = defaultdict(list)  # {branch_id : [branch_segment_ids]}
    branch_id_to_vessel_id_map = defaultdict(
        dict
    )  # {branch_id : {branch_segment_id : vessel_id}}
    vessel_id_to_length_map = get_vessel_id_to_length_map(parameters)

    for vessel in parameters["vessels"]:
        vessel_id = vessel["vessel_id"]
        vessel_name = vessel["vessel_name"]
        vessel_name_split = vessel_name.split("_")
        branch_id = int(
            (re.match(r"([a-z]+)([0-9]+)", vessel_name_split[0], re.I)).groups()[1]
        )
        branch_segment_id = int(
            (re.match(r"([a-z]+)([0-9]+)", vessel_name_split[1], re.I)).groups()[1]
        )
        branch_segment_ids[branch_id].append(branch_segment_id)
        branch_id_to_vessel_id_map[branch_id][branch_segment_id] = vessel_id

    for branch_id in branch_segment_ids:
        num_nodes_for_branch = max(branch_segment_ids[branch_id]) + 2
        for qoi in zero_d_results:
            zero_d_results[qoi][branch_id] = np.zeros(
                (num_nodes_for_branch, num_time_pts)
            )
        zero_d_results["distance"][branch_id] = np.zeros(num_nodes_for_branch)
        branch_segment_ids[
            branch_id
        ].sort()  # sort by ascending order of branch_segment_id
        counter = 0
        for branch_segment_id in branch_segment_ids[branch_id]:
            vessel_id = branch_id_to_vessel_id_map[branch_id][branch_segment_id]
            zero_d_results["distance"][branch_id][branch_segment_id + 1] = (
                zero_d_results["distance"][branch_id][counter]
                + vessel_id_to_length_map[vessel_id]
            )
            counter += 1
    zero_d_results["time"] = zero_d_time
    return zero_d_results


def get_vessel_id_to_vessel_name_map(parameters):
    vessel_id_to_vessel_name_map = {}
    for vessel in parameters["vessels"]:
        vessel_id = vessel["vessel_id"]
        vessel_name = vessel["vessel_name"]
        vessel_id_to_vessel_name_map[vessel_id] = vessel_name
    return vessel_id_to_vessel_name_map


def reformat_network_util_results_branch(
    zero_d_time, results_0d, var_name_list, parameters
):
    """
    Purpose:
        Reformat the 0d simulation results for just the branches into a dictionary (zero_d_results)
    Inputs:
        np.array zero_d_time
            = np.array of simulated time points
        np.array results_0d
            = np.array of the 0d simulation results, where the rows correspond to the each simulated time point and column j corresponds to the 0d solution for the solution variable name in var_name_list[j]
        list var_name_list
            = list of the 0d simulation results' solution variable names; most of the items in var_name_list are the QoIs + the names of the wires used in the 0d model (the wires connecting the 0d LPNBlock objects), where the wire names are usually comprised of the wire's inlet block name + "_" + the wire's outlet block name
                Example:
                    for var_name_list = ['P_V6_BC6_outlet', 'Q_V6_BC6_outlet'], then results_0d[:, i] holds the pressure (i = 0) or flow rate simulation result (i = 1) (both as np.arrays) for wire R6_BC6_outlet. This wire connects a resistance vessel block to an outlet BC block (specifically for vessel segment #6)
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
    Returns:
        dict zero_d_results
            =   {
                    "time" : 1d np.array of simulated time points,

                    "distance" : {branch_id : 1d np.array of distance of the branch's 0d nodes along the centerline}

                    "flow" : {branch_id : 2d np.array of flow rate where each row represents a 0d node on the branch and each column represents a time point,

                    "pressure" : {branch_id : 2d np.array of pressure where each row represents a 0d node on the branch and each column represents a time point}
                }

            - examples:
                1. plt.plot(zero_d_results["time"], zero_d_results["flow"][branch_id][0, :])
                        --> plot time vs the 0d flow waveform for the 0th node on branch, branch_id; this yields a plot that shows how the flow rate changes over time

                2. plt.plot(zero_d_results["distance"], zero_d_results["pressure"][branch_id][:, -1])
                        --> plot centerline distance vs the 0d pressure (at the last simulated time step); this yields a plot that shows how the pressure changes along the axial dimension of a vessel
    """
    qoi_map = {"Q": "flow", "P": "pressure"}
    zero_d_results = initialize_0d_results_dict_branch(parameters, zero_d_time)
    vessel_id_to_vessel_name_map = get_vessel_id_to_vessel_name_map(parameters)

    for i in range(len(var_name_list)):
        if "var" not in var_name_list[i]:  # var_name_list[i] == wire_name
            # the only possible combination of wire connections are: 1) vessel <--> vessel 2) vessel <--> junction 3) vessel <--> boundary condition; in all of these cases, there is at least one vessel block ("V") in each wire_name
            if "V" not in var_name_list[i]:
                message = "Error. It is expected that every wire in the 0d model must be connected to at least one vessel block."
                raise RuntimeError(message)
            else:
                # parse var_name
                var_name_split = var_name_list[i].split("_")
                qoi_header = var_name_split[0]

                if var_name_split[1].startswith(
                    "V"
                ):  # the wire connected downstream of this vessel block
                    vessel_id = int(var_name_split[1][1:])
                    vessel_name = vessel_id_to_vessel_name_map[vessel_id]
                    vessel_name_split = vessel_name.split("_")
                    branch_id = int(
                        (
                            re.match(r"([a-z]+)([0-9]+)", vessel_name_split[0], re.I)
                        ).groups()[1]
                    )
                    branch_segment_id = int(
                        (
                            re.match(r"([a-z]+)([0-9]+)", vessel_name_split[1], re.I)
                        ).groups()[1]
                    )
                    branch_node_id = branch_segment_id + 1
                    zero_d_results[qoi_map[qoi_header]][branch_id][
                        branch_node_id, :
                    ] = results_0d[:, i]
                else:  # need to find the inlet wire/node of the branch
                    if var_name_split[1].startswith(
                        "BC"
                    ):  # inlet wire/node of the branch
                        vessel_id = int(var_name_split[3][1:])
                        vessel_name = vessel_id_to_vessel_name_map[vessel_id]
                        vessel_name_split = vessel_name.split("_")
                        branch_id = int(
                            (
                                re.match(
                                    r"([a-z]+)([0-9]+)", vessel_name_split[0], re.I
                                )
                            ).groups()[1]
                        )
                        branch_segment_id = int(
                            (
                                re.match(
                                    r"([a-z]+)([0-9]+)", vessel_name_split[1], re.I
                                )
                            ).groups()[1]
                        )
                        if branch_segment_id != 0:
                            message = "Error. branch_segment_id should be 0 here because we are at the inlet wire of the branch."
                            raise RuntimeError(message)
                        else:
                            branch_node_id = branch_segment_id
                            zero_d_results[qoi_map[qoi_header]][branch_id][
                                branch_node_id, :
                            ] = results_0d[:, i]
                    elif var_name_split[1].startswith(
                        "J"
                    ):  # this wire could either be 1) the inlet wire/node of a branch or 2) some internal wire in the branch (where that internal wire is 1) connecting 2 vessel blocks or 2) connecting a vessel block and a junction block) (and we dont care about internal wires)
                        vessel_id = int(var_name_split[2][1:])
                        vessel_name = vessel_id_to_vessel_name_map[vessel_id]
                        vessel_name_split = vessel_name.split("_")
                        branch_id = int(
                            (
                                re.match(
                                    r"([a-z]+)([0-9]+)", vessel_name_split[0], re.I
                                )
                            ).groups()[1]
                        )
                        branch_segment_id = int(
                            (
                                re.match(
                                    r"([a-z]+)([0-9]+)", vessel_name_split[1], re.I
                                )
                            ).groups()[1]
                        )
                        if (
                            branch_segment_id == 0
                        ):  # this is the inlet wire/node of the branch
                            branch_node_id = branch_segment_id
                            zero_d_results[qoi_map[qoi_header]][branch_id][
                                branch_node_id, :
                            ] = results_0d[:, i]
                        else:  # this is an internal wire/node in the branch
                            pass  # do nothing here, since we are ignoring the internal wires where the vessel block is connected downstream of the wire (and the junction block is connected upstream of the wire)
                    else:
                        message = 'Error. It is not possible for a block name to begin with something other than, "V", "J", or "BC".'
                        raise RuntimeError(message)

    return zero_d_results


def extract_last_cardiac_cycle_simulation_results(
    time, results, number_of_time_pts_per_cardiac_cycle
):
    """
    Purpose:
        Extract the simulation results for the last cardiac cycle for the given np.arrays, time and results
    Inputs:
        np.array time
            = array of time points
        np.array results
            = 2d array of simulation results
        int number_of_time_pts_per_cardiac_cycle
            = number of simulated time points per cycle
    Returns:
        np.array time_for_last_cardiac_cycle
            = time, but condensed to contain just the values for the last cardiac cycle
        np.array results_for_last_cardiac_cycle
            = results, but condensed to contain just the values for the last cardiac cycle
    """
    time_for_last_cardiac_cycle = time[-1 * number_of_time_pts_per_cardiac_cycle :]
    time_for_last_cardiac_cycle = (
        time_for_last_cardiac_cycle - time_for_last_cardiac_cycle[0] + time[0]
    )
    return (
        time_for_last_cardiac_cycle,
        results[-1 * number_of_time_pts_per_cardiac_cycle :, :],
    )


def run_last_cycle_extraction_routines(
    cardiac_cycle_period, number_of_time_pts_per_cardiac_cycle, zero_d_time, results_0d
):
    """
    Purpose:
        Extract the last cardiac cycle waveform for all 0d simulation results
    Inputs:
        float cardiac_cycle_period
            = period of a cardiac cycle
        int number_of_time_pts_per_cardiac_cycle
            = number of simulated 0d time points per cycle
        np.array zero_d_time
            = np.array of simulated time points
        np.array results_0d
            = np.array of the 0d simulation results, where the rows correspond to the each simulated time point and column j corresponds to the 0d solution for the solution variable name in var_name_list[j]
    Returns:
        np.array time_for_last_cardiac_cycle
            = time, but condensed to contain just the values for the last cardiac cycle
        np.array res
            = results, but condensed to contain just the values for the last cardiac cycle
    """

    time_for_last_cardiac_cycle, res = extract_last_cardiac_cycle_simulation_results(
        zero_d_time, results_0d, number_of_time_pts_per_cardiac_cycle
    )

    # check that the cardiac cycle period is correctly given by zero_d_time
    errr = (
        (
            cardiac_cycle_period
            - (time_for_last_cardiac_cycle[-1] - time_for_last_cardiac_cycle[0])
        )
        / cardiac_cycle_period
        * 100.0
    )
    if errr > 0.01:
        message = "Error. cardiac_cycle_period != time_for_last_cardiac_cycle[-1] - time_for_last_cardiac_cycle[0]"
        raise RuntimeError(message)

    return time_for_last_cardiac_cycle, res


def save_simulation_results(zero_d_simulation_results_file_path, zero_d_results):
    """
    Purpose:
        Save the 0d simulation results to a .npy file.
    Usage:
        To open and load the .npy file to extract the 0d simulation results, use the following command:
            zero_d_results = np.load(/path/to/zero/d/simulation/results/npy/file, allow_pickle = True).item()
    Inputs:
        string zero_d_simulation_results_file_path
            = path to the .npy file to which the 0d simulation results will be saved
        dict zero_d_results
            = obtained from reformat_network_util_results_all or reformat_network_util_results_branch
    Returns:
        void, but saves a .npy file storing the 0d simulation results (zero_d_results) as a dictionary.
    """
    np.save(zero_d_simulation_results_file_path, zero_d_results)


def set_up_and_run_0d_simulation(
    zero_d_solver_input_file_path,
    last_cycle=False,
    save_results_all=False,
    save_results_branch=True,
    use_ICs_from_npy_file=False,
    ICs_npy_file_path=None,
    save_y_ydot_to_npy=False,
    y_ydot_file_path=None,
    simulation_start_time=0.0,
    use_steady_soltns_as_ics=True,
):
    """
    Purpose:
        Create all network_util_NR::LPNBlock objects for the 0d model and run the 0d simulation.
    Inputs:
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
        boolean draw_directed_graph
            = True to visualize the 0d model as a directed graph using networkx -- saves the graph to a .png file (hierarchical graph layout) and a networkx .dot file; False, otherwise. .dot file can be opened with neato from graphviz to visualize the directed in a different format.
        boolean last_cycle
            = True to return 0d simulation results for only the last cycle; False to return the results for all simulated cycles
        boolean save_results_all
            = True to save all 0d simulation results (reformatted to a readable dictionary format) to a .npy file
        boolean save_results_branch
            = True to save the 0d simulation results for just the branches (reformatted to a readable dictionary format, while preserving the 1d/centerline branch structure) to a .npy file
        boolean use_custom_0d_elements
            = True to use user-defined, custom 0d elements in the 0d model; False, otherwire
        boolean use_ICs_from_npy_file
            = True to use user-prescribed ICs (saved in a .npy file)
        string ICs_npy_file_path
            = path to .npy file storing the initial conditions
        boolean save_y_ydot_to_npy
            = True to save the entire 0d solution and its time derivative at the final time step to a .npy file
        string y_ydot_file_path
            = name of the .npy file to which the 0d solution and its time derivative at the final time step will be saved
        float simulation_start_time
            = time at which to begin the 0d simulation
            -- assumes the cardiac cycle begins at t = 0.0 and ends at t = cardiac_cycle_period
            -- 0.0 <= simulation_start_time <= cardiac_cycle_period
        boolean use_steady_soltns_as_ics
            = True to 1) run the 0d simulation using steady mean BCs and zero initial conditions and then 2) use the steady-state solution from the previous part as the initial condition for the original pulsatile simulation
    Caveats:
        The save_results_branch option works only for 0d models with the branching structure where each vessel is modeled as a single branch with 1 or multiple sub-segments
    Returns:
        void
    """

    with open(zero_d_solver_input_file_path, "r") as infile:
        parameters = json.load(infile)

    if use_steady_soltns_as_ics:
        parameters_mean = copy.deepcopy(parameters)
        (
            parameters_mean,
            altered_bc_blocks,
        ) = use_steady_bcs.convert_unsteady_bcs_to_steady(parameters_mean)

        # to run the 0d model with steady BCs to steady-state, simulate this model with large time step size for an arbitrarily small number of cardiac cycles
        parameters_mean["simulation_parameters"][
            "number_of_time_pts_per_cardiac_cycle"
        ] = 11
        parameters_mean["simulation_parameters"]["number_of_cardiac_cycles"] = 3

        y_ydot_file_path_temp = (
            os.path.splitext(zero_d_solver_input_file_path)[0]
            + "_initial_conditions.npy"
        )

        create_LPN_blocks(parameters_mean)
        set_solver_parameters(parameters_mean)
        (
            zero_d_time,
            results_0d,
            var_name_list,
            y_f,
            ydot_f,
            var_name_list_original,
        ) = run_network_util(
            parameters_mean,
            use_ICs_from_npy_file=False,
            ICs_npy_file_path=None,
            save_y_ydot_to_npy=False,
            y_ydot_file_path=None,
            simulation_start_time=simulation_start_time,
        )

        (
            y0,
            ydot0,
            var_name_list,
        ) = use_steady_bcs.restore_internal_variables_for_capacitance_based_bcs(
            y_f, ydot_f, var_name_list_original, altered_bc_blocks
        )

        save_ics(y_ydot_file_path_temp, y0, ydot0, var_name_list)

        use_ICs_from_npy_file = True
        ICs_npy_file_path = y_ydot_file_path_temp

    create_LPN_blocks(parameters)
    set_solver_parameters(parameters)
    zero_d_time, results_0d, var_name_list, _, _, _ = run_network_util(
        parameters,
        use_ICs_from_npy_file,
        ICs_npy_file_path,
        save_y_ydot_to_npy,
        y_ydot_file_path,
        simulation_start_time,
    )
    print("0D simulation completed!\n")

    # postprocessing
    if last_cycle == True:
        zero_d_time, results_0d = run_last_cycle_extraction_routines(
            parameters["simulation_parameters"]["cardiac_cycle_period"],
            parameters["simulation_parameters"]["number_of_time_pts_per_cardiac_cycle"],
            zero_d_time,
            results_0d,
        )
    if save_results_all or save_results_branch:
        zero_d_input_file_name = os.path.splitext(zero_d_solver_input_file_path)[0]
        if save_results_all:
            zero_d_simulation_results_file_path = (
                zero_d_input_file_name + "_all_results"
            )
            zero_d_results = reformat_network_util_results_all(
                zero_d_time, results_0d, var_name_list
            )
            save_simulation_results(zero_d_simulation_results_file_path, zero_d_results)
        if save_results_branch:
            zero_d_simulation_results_file_path = (
                zero_d_input_file_name + "_branch_results"
            )
            zero_d_results = reformat_network_util_results_branch(
                zero_d_time, results_0d, var_name_list, parameters
            )
            save_simulation_results(zero_d_simulation_results_file_path, zero_d_results)


def create_parser():
    """Create a command-line parser."""

    # Create the parser.
    parser = argparse.ArgumentParser(description="This code runs the 0d solver.")

    # Define 0D solver commands.

    parser.add_argument("zero", help="Path to 0d solver input file")

    parser.add_argument(
        "-v",
        "--visualize",
        action="store_true",
        help="Visualize the 0d model as a networkx directed graph and save to .png file",
    )

    parser.add_argument(
        "-l",
        "--returnLast",
        action="store_true",
        help="Return results for only the last simulated cardiac cycle",
    )

    parser.add_argument(
        "-sa",
        "--saveAll",
        default=True,
        action="store_true",
        help="Save all simulation results to a .npy file",
    )

    parser.add_argument(
        "-sb",
        "--saveBranch",
        default=True,
        action="store_true",
        help="Save the simulation results (preserving the 1d/centerline branch structure) to a .npy file",
    )  # todo: do we need to change action to 'store_false' here?

    parser.add_argument(
        "-i",
        "--useICs",
        action="store_true",
        help="Use initial conditions from .npy file",
    )  # todo: need to prevent users from using both: useSteadyIC and useICs

    parser.add_argument(
        "-pi",
        "--ICsPath",
        help="Path to the .npy file containing the initial conditions",
    )

    parser.add_argument(
        "-y", "--saveYydot", action="store_true", help="Save y and ydot to a .npy file"
    )

    parser.add_argument(
        "-py", "--yydotPath", help="Path to the .npy file containing y and ydot"
    )

    parser.add_argument(
        "-j", "--checkJacobian", action="store_true", help="Check the Jacobian"
    )

    parser.add_argument(
        "-it",
        "--initialTime",
        default=0.0,
        type=float,
        help="Start (initial) time of the 0d simulation",
    )

    parser.add_argument(
        "-sic",
        "--useSteadyIC",
        action="store_true",
        help="Run the pulsatile 0d simulation using the steady-state solution from the equivalent steady 0d model as the initial conditions.",
    )  # caveat - does not work with custom, user-defined BCs

    return parser


def main():

    args = sys.argv[1:]

    parser = create_parser()

    args = parser.parse_args(args)

    set_up_and_run_0d_simulation(
        zero_d_solver_input_file_path=args.zero,
        last_cycle=args.returnLast,
        save_results_all=args.saveAll,
        save_results_branch=args.saveBranch,
        use_custom_0d_elements=args.useCustom,
        use_ICs_from_npy_file=args.useICs,
        ICs_npy_file_path=args.ICsPath,
        save_y_ydot_to_npy=args.saveYydot,
        y_ydot_file_path=args.yydotPath,
        simulation_start_time=args.initialTime,
        use_steady_soltns_as_ics=args.useSteadyIC,
    )


if __name__ == "__main__":
    main(main)
