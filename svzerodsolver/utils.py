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

from copy import deepcopy

import numpy as np

from .model.bloodvessel import BloodVessel
from .model.dofhandler import DOFHandler
from .model.flowreferencebc import FlowReferenceBC
from .model.junction import Junction
from .model.node import Node
from .model.openloopcoronarybc import OpenLoopCoronaryBC
from .model.pressurereferencebc import PressureReferenceBC
from .model.resistancebc import ResistanceBC
from .model.windkesselbc import WindkesselBC


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


def create_vessel_blocks(parameters, steady):
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
        if "boundary_conditions" in vessel_config:
            locations = [
                loc
                for loc in ("inlet", "outlet")
                if loc in vessel_config["boundary_conditions"]
            ]
            for location in locations:
                vessel_id = vessel_config["vessel_id"]
                for config in parameters["boundary_conditions"]:
                    if (
                        config["bc_name"]
                        == vessel_config["boundary_conditions"][location]
                    ):
                        bc_config = dict(
                            name="BC" + str(vessel_id) + "_" + location,
                            steady=steady,
                            **config,
                        )
                        break

                if (
                    "t" in bc_config["bc_values"]
                    and len(bc_config["bc_values"]["t"]) >= 2
                ):
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
                    elif (
                        "cardiac_cycle_period"
                        not in parameters["simulation_parameters"]
                    ):
                        parameters["simulation_parameters"][
                            "cardiac_cycle_period"
                        ] = cardiac_cycle_period

                if bc_config["bc_type"] == "RESISTANCE":
                    bc = ResistanceBC.from_config(bc_config)
                elif bc_config["bc_type"] == "RCR":
                    bc = WindkesselBC.from_config(bc_config)
                elif bc_config["bc_type"] == "FLOW":
                    bc = FlowReferenceBC.from_config(bc_config)
                elif bc_config["bc_type"] == "PRESSURE":
                    bc = PressureReferenceBC.from_config(bc_config)
                elif bc_config["bc_type"] == "CORONARY":
                    bc = OpenLoopCoronaryBC.from_config(bc_config)
                else:
                    raise NotImplementedError
                block_dict[bc.name] = bc
    return block_dict, connections


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
    vessel_blocks, vessel_connections = create_vessel_blocks(parameters, steady=steady)
    all_blocks = junction_blocks | vessel_blocks
    dofhandler = DOFHandler()
    for ele1_name, ele2_name in junction_connections + vessel_connections:
        node = Node(
            all_blocks[ele1_name],
            all_blocks[ele2_name],
            name=ele1_name + "_" + ele2_name,
        )
        node.setup_dofs(dofhandler)
    for key in sorted(all_blocks):
        all_blocks[key].setup_dofs(dofhandler)
    return list(all_blocks.values()), dofhandler


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
