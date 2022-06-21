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
from . import use_steady_bcs
from .results import reformat_network_util_results_branch


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

        block_list, dofhandler = create_blocks(parameters_mean)
        set_solver_parameters(parameters_mean)
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
            y_list[-1],
            ydot_list[-1],
            list(dofhandler.variables.values()),
            altered_bc_blocks,
        )

    block_list, dofhandler = create_blocks(parameters)
    y_initial, ydot_initial = convert_ics(
        list(dofhandler.variables.values()),
        var_name_list_original,
        y_initial,
        ydot_initial,
    )
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
