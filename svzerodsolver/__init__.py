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
import numpy as np

from .blocks import create_blocks, BloodVessel
from .integrator import run_integrator

import numpy as np
from copy import deepcopy

NAME = "svZeroDSolver"
VERSION = "v1.0"
COPYRIGHT = (
    "Stanford University, The Regents of the University of California, and others."
)


def convert_unsteady_bcs_to_steady(parameters):
    steady_parameters = deepcopy(parameters)
    steady_parameters["simulation_parameters"][
        "number_of_time_pts_per_cardiac_cycle"
    ] = 11
    steady_parameters["simulation_parameters"]["number_of_cardiac_cycles"] = 3
    bc_identifiers = {"FLOW": "Q", "PRESSURE": "P", "CORONARY": "Pim"}

    for i, bc in enumerate(parameters["boundary_conditions"]):
        if bc["bc_type"] in bc_identifiers:
            bc_values = bc["bc_values"][bc_identifiers[bc["bc_type"]]]
            # Time averaged value for a single cariadic_cycle
            del steady_parameters["boundary_conditions"][i]["bc_values"]["t"]
            steady_parameters["boundary_conditions"][i]["bc_values"][
                bc_identifiers[bc["bc_type"]]
            ] = np.mean(bc_values)
        if bc["bc_type"] == "RCR":
            steady_parameters["boundary_conditions"][i]["bc_values"]["C"] = 0.0

    return steady_parameters


def get_solver_params(parameters):
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
    sim_params = parameters["simulation_parameters"]
    cardiac_cycle_period = sim_params.get("cardiac_cycle_period", 1.0)
    num_cycles = sim_params.get("number_of_cardiac_cycles")
    num_pts_per_cycle = sim_params.get("number_of_time_pts_per_cardiac_cycle")
    time_step_size = cardiac_cycle_period / (num_pts_per_cycle - 1)
    num_time_steps = int((num_pts_per_cycle - 1) * num_cycles + 1)
    return time_step_size, num_time_steps


def format_results_to_dict(zero_d_time, results_0d, block_list):

    results_0d = np.array(results_0d)

    vessels = [block for block in block_list if isinstance(block, BloodVessel)]

    results = {
        "flow_in": [],
        "flow_out": [],
        "names": [],
        "pressure_in": [],
        "pressure_out": [],
        "time": list(zero_d_time),
    }

    for vessel in vessels:

        results["names"].append(vessel.name)
        results["flow_in"].append(list(results_0d[:, vessel.inflow_nodes[0].flow_dof]))
        results["flow_out"].append(
            list(results_0d[:, vessel.outflow_nodes[0].flow_dof])
        )
        results["pressure_in"].append(
            list(results_0d[:, vessel.inflow_nodes[0].pres_dof])
        )
        results["pressure_out"].append(
            list(results_0d[:, vessel.outflow_nodes[0].pres_dof])
        )

    return results


def run_simulation_from_config(
    parameters,
    use_steady_soltns_as_ics=True,
):

    y_initial = None
    ydot_initial = None
    if use_steady_soltns_as_ics:
        steady_parameters = convert_unsteady_bcs_to_steady(parameters)

        # to run the 0d model with steady BCs to steady-state, simulate this model with large time step size for an arbitrarily small number of cardiac cycles
        block_list, dofhandler = create_blocks(steady_parameters, steady=True)
        time_step_size, num_time_steps = get_solver_params(steady_parameters)
        (_, y_list, ydot_list) = run_integrator(
            block_list,
            dofhandler,
            num_time_steps,
            time_step_size,
        )
        y_initial = y_list[-1]
        ydot_initial = ydot_list[-1]

    block_list, dofhandler = create_blocks(parameters)
    time_step_size, num_time_steps = get_solver_params(parameters)
    (time_steps, y_list, ydot_list) = run_integrator(
        block_list,
        dofhandler,
        num_time_steps,
        time_step_size,
        y_initial,
        ydot_initial,
    )

    zero_d_results_branch = format_results_to_dict(time_steps, y_list, block_list)
    return zero_d_results_branch
