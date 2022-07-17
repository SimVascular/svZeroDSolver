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
"""This module holds helper functions for svZeroDSolver."""
from copy import deepcopy

import numpy as np

from .model.bloodvessel import BloodVessel
from .model.dofhandler import DOFHandler
from .model.flowreferencebc import FlowReferenceBC
from .model.internaljunction import InternalJunction
from .model.node import Node
from .model.openloopcoronarybc import OpenLoopCoronaryBC
from .model.pressurereferencebc import PressureReferenceBC
from .model.resistancebc import ResistanceBC
from .model.windkesselbc import WindkesselBC


def get_solver_params(config):
    """Determine time step size and number of timesteps.

    Args:
        config: svZeroDSolver configuration.

    Returns:
        time_step_size: Time step size:
        num_time_steps: Number of time steps.
    """
    sim_params = config["simulation_parameters"]
    cardiac_cycle_period = sim_params.get("cardiac_cycle_period", 1.0)
    num_cycles = sim_params.get("number_of_cardiac_cycles")
    num_pts_per_cycle = sim_params.get("number_of_time_pts_per_cardiac_cycle")
    time_step_size = cardiac_cycle_period / (num_pts_per_cycle - 1)
    num_time_steps = int((num_pts_per_cycle - 1) * num_cycles + 1)
    return time_step_size, num_time_steps


def convert_unsteady_bcs_to_steady(config):
    """Convert unsteady boundary conditions to their steady equivalent.

    Args:
        config: svZeroDSolver configuration.

    Returns:
        steady_config: Configuration of the steady equivalent.
    """
    steady_config = deepcopy(config)
    steady_config["simulation_parameters"][
        "number_of_time_pts_per_cardiac_cycle"
    ] = 11
    steady_config["simulation_parameters"]["number_of_cardiac_cycles"] = 3
    bc_identifiers = {"FLOW": "Q", "PRESSURE": "P", "CORONARY": "Pim"}
    for i, bc in enumerate(config["boundary_conditions"]):
        if bc["bc_type"] in bc_identifiers:
            bc_values = bc["bc_values"][bc_identifiers[bc["bc_type"]]]
            # Time averaged value for a single cariadic_cycle
            del steady_config["boundary_conditions"][i]["bc_values"]["t"]
            steady_config["boundary_conditions"][i]["bc_values"][
                bc_identifiers[bc["bc_type"]]
            ] = np.mean(bc_values)
        if bc["bc_type"] == "RCR":
            steady_config["boundary_conditions"][i]["bc_values"]["C"] = 0.0

    return steady_config


def create_blocks(config, steady=False):
    """Create blocks.

    Args:
        config: svZeroDSolver configuration.
        steady: Toggle if blocks should be created with steady behavior.

    Returns:
        blocks: List of created blocks.
        dofhandler: Degree-of-freedom handler.
    """
    block_dict = {}
    connections = []

    # Create junctions
    for junction_config in config["junctions"]:
        if junction_config["junction_type"] in [
            "NORMAL_JUNCTION",
            "internal_junction",
        ]:
            junction = InternalJunction.from_config(junction_config)
        else:
            raise ValueError(
                f"Unknown junction type: {junction_config['junction_type']}"
            )
        connections += [
            (f"V{vid}", junction.name)
            for vid in junction_config["inlet_vessels"]
        ]
        connections += [
            (junction.name, f"V{vid}")
            for vid in junction_config["outlet_vessels"]
        ]
        if junction.name in block_dict:
            raise RuntimeError(f"Junction {junction.name} already exists.")
        block_dict[junction.name] = junction

    # Create vessels and boundary conditions
    for vessel_config in config["vessels"]:
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
                for bc_config in config["boundary_conditions"]:
                    if (
                        bc_config["bc_name"]
                        == vessel_config["boundary_conditions"][location]
                    ):
                        bc_config = dict(
                            name="BC" + str(vessel_id) + "_" + location,
                            steady=steady,
                            **bc_config,
                        )
                        break

                if (
                    "t" in bc_config["bc_values"]
                    and len(bc_config["bc_values"]["t"]) >= 2
                ):
                    cardiac_cycle_period = (
                        bc_config["bc_values"]["t"][-1]
                        - bc_config["bc_values"]["t"][0]
                    )
                    if (
                        "cardiac_cycle_period"
                        in config["simulation_parameters"]
                        and cardiac_cycle_period
                        != config["simulation_parameters"][
                            "cardiac_cycle_period"
                        ]
                    ):
                        raise RuntimeError(
                            f"The time series of the boundary condition for segment {vessel_id} does "
                            "not have the same cardiac cycle period as the other boundary conditions."
                        )
                    elif (
                        "cardiac_cycle_period"
                        not in config["simulation_parameters"]
                    ):
                        config["simulation_parameters"][
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

    # Create nodes
    dofhandler = DOFHandler()
    for ele1_name, ele2_name in connections:
        node = Node(
            block_dict[ele1_name],
            block_dict[ele2_name],
            name=ele1_name + "_" + ele2_name,
        )
        node.setup_dofs(dofhandler)

    # Setup degrees of freedom of the model
    for key in sorted(block_dict):
        block_dict[key].setup_dofs(dofhandler)
    return list(block_dict.values()), dofhandler


def format_results_to_dict(time_steps, result_array, block_list):
    """Format result array to dict format.

    Args:
        time_steps: List of time steps.

    Returns:
        results: Resuluts in dict format.
    """

    result_array = np.array(result_array)
    vessels = [block for block in block_list if isinstance(block, BloodVessel)]
    results = {
        "flow_in": [],
        "flow_out": [],
        "names": [],
        "pressure_in": [],
        "pressure_out": [],
        "time": list(time_steps),
    }

    for vessel in vessels:

        results["names"].append(vessel.name)
        results["flow_in"].append(
            list(result_array[:, vessel.inflow_nodes[0].flow_dof])
        )
        results["flow_out"].append(
            list(result_array[:, vessel.outflow_nodes[0].flow_dof])
        )
        results["pressure_in"].append(
            list(result_array[:, vessel.inflow_nodes[0].pres_dof])
        )
        results["pressure_out"].append(
            list(result_array[:, vessel.outflow_nodes[0].pres_dof])
        )

    return results
