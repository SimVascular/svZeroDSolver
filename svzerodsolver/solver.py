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
import copy
import numpy as np

from .blocks import create_blocks
from .time_integration import run_integrator
from .results import reformat_network_util_results_branch

import numpy as np
from copy import deepcopy


def convert_unsteady_bcs_to_steady(parameters):
    steady_parameters = deepcopy(parameters)
    bc_identifiers = {"FLOW": "Q", "PRESSURE": "P", "CORONARY": "Pim"}

    for i, bc in enumerate(parameters["boundary_conditions"]):
        if bc["bc_type"] in bc_identifiers:
            time = bc["bc_values"]["t"]
            bc_values = bc["bc_values"][bc_identifiers[bc["bc_type"]]]
            # Time averaged value for a single cariadic_cycle
            time_averaged_value = 1.0 / (time[-1] - time[0]) * np.trapz(bc_values, time)
            steady_parameters["boundary_conditions"][i]["bc_values"]["t"] = [
                time[0],
                time[-1],
            ]
            steady_parameters["boundary_conditions"][i]["bc_values"][
                bc_identifiers[bc["bc_type"]]
            ] = [time_averaged_value, time_averaged_value]

    return steady_parameters


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


def run_simulation_from_config(
    parameters,
    use_steady_soltns_as_ics=True,
):

    y_initial = None
    ydot_initial = None
    if use_steady_soltns_as_ics:
        steady_paramters = convert_unsteady_bcs_to_steady(parameters)

        # to run the 0d model with steady BCs to steady-state, simulate this model with large time step size for an arbitrarily small number of cardiac cycles
        steady_paramters["simulation_parameters"][
            "number_of_time_pts_per_cardiac_cycle"
        ] = 11
        steady_paramters["simulation_parameters"]["number_of_cardiac_cycles"] = 3

        block_list, dofhandler = create_blocks(steady_paramters, steady=True)
        set_solver_parameters(steady_paramters)
        (_, y_list, ydot_list) = run_integrator(
            block_list,
            steady_paramters["simulation_parameters"][
                "total_number_of_simulated_time_steps"
            ],
            steady_paramters["simulation_parameters"]["delta_t"],
        )
        y_initial = y_list[-1]
        ydot_initial = ydot_list[-1]

    block_list, dofhandler = create_blocks(parameters)
    set_solver_parameters(parameters)
    (time_steps, y_list, ydot_list) = run_integrator(
        block_list,
        parameters["simulation_parameters"]["total_number_of_simulated_time_steps"],
        parameters["simulation_parameters"]["delta_t"],
        y_initial,
        ydot_initial,
    )

    zero_d_results_branch = reformat_network_util_results_branch(
        time_steps, np.array(y_list), list(dofhandler.variables.values()), parameters
    )
    return zero_d_results_branch
