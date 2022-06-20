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
import re
import copy
import numpy as np

from .blocks import create_LPN_blocks
from . import connections
from .time_integration import run_integrator
from . import use_steady_bcs

from collections import defaultdict


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


def convert_ics(var_name_list, var_name_list_original, y, ydot):

    var_name_list_loaded = var_name_list_original
    y_initial_loaded = y
    ydot_initial_loaded = ydot

    y_initial = np.zeros(len(var_name_list))
    ydot_initial = np.zeros(len(var_name_list))
    for i in range(len(var_name_list)):
        var_name = var_name_list[i]
        ind = var_name_list_loaded.index(var_name)
        y_initial[i] = y_initial_loaded[ind]
        ydot_initial[i] = ydot_initial_loaded[ind]

    return y_initial, ydot_initial


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


def run_simulation_from_config(
    parameters,
    use_steady_soltns_as_ics=True,
):

    y_initial = None
    ydot_initial = None
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

        block_list = create_LPN_blocks(parameters_mean)
        set_solver_parameters(parameters_mean)
        wire_dict = connections.connect_blocks_by_inblock_list(block_list)
        var_name_list_original = connections.assign_global_ids(block_list, wire_dict)
        (_, y_list, ydot_list) = run_integrator(
            block_list,
            parameters_mean["simulation_parameters"][
                "total_number_of_simulated_time_steps"
            ],
            parameters_mean["simulation_parameters"]["delta_t"],
        )

        (
            y_initial,
            ydot_initial,
            var_name_list_original,
        ) = use_steady_bcs.restore_internal_variables_for_capacitance_based_bcs(
            y_list[-1], ydot_list[-1], var_name_list_original, altered_bc_blocks
        )

    block_list = create_LPN_blocks(parameters)
    set_solver_parameters(parameters)
    wire_dict = connections.connect_blocks_by_inblock_list(block_list)
    var_name_list = connections.assign_global_ids(block_list, wire_dict)
    y_initial, ydot_initial = convert_ics(
        var_name_list, var_name_list_original, y_initial, ydot_initial
    )
    (time_steps, y_list, ydot_list) = run_integrator(
        block_list,
        parameters["simulation_parameters"]["total_number_of_simulated_time_steps"],
        parameters["simulation_parameters"]["delta_t"],
        y_initial,
        ydot_initial,
    )

    zero_d_results_all = reformat_network_util_results_all(
        time_steps, np.array(y_list), var_name_list
    )

    zero_d_results_branch = reformat_network_util_results_branch(
        time_steps, np.array(y_list), var_name_list, parameters
    )
    return zero_d_results_branch, zero_d_results_all
