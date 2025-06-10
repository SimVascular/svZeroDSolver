# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
import copy
import numpy as np

def get_bc_name_to_index_map(parameters):
    bc_name_to_index_map = {}
    for i in range(len(parameters["boundary_conditions"])):
        bc_name = parameters["boundary_conditions"][i]["bc_name"]
        bc_name_to_index_map[bc_name] = i
    return bc_name_to_index_map


def create_block_to_boundary_condition_map(parameters):
    def add_bc_to_map(block_name, bc_name, location):
        if block_name not in block_to_boundary_condition_map:
            block_to_boundary_condition_map[block_name] = {}
        for boundary_condition in parameters["boundary_conditions"]:
            if boundary_condition["bc_name"] == bc_name:
                block_to_boundary_condition_map[block_name][location] = boundary_condition

    block_to_boundary_condition_map = {}

    # Process vessels
    for vessel in parameters["vessels"]:
        if "boundary_conditions" in vessel:
            vessel_id = vessel["vessel_id"]
            for location, bc_name in vessel["boundary_conditions"].items():
                add_bc_to_map(vessel_id, bc_name, location)

    # Process valves
    if 'valves' in parameters:
        for valve in parameters['valves']:
            valve_name = valve['name']
            upstream_bc = valve['params']['upstream_block']
            downstream_bc = valve['params']['downstream_block']

            if upstream_bc in get_bc_name_to_index_map(parameters):
                add_bc_to_map(valve_name, upstream_bc, 'inlet')

            if downstream_bc in get_bc_name_to_index_map(parameters):
                add_bc_to_map(valve_name, downstream_bc, 'outlet')

    parameters["block_to_boundary_condition_map"] = block_to_boundary_condition_map



def get_ids_of_cap_vessels(parameters, location):  # location == "inlet" or "outlet"
    ids_of_cap_vessels = []
    for vessel in parameters["vessels"]:
        if "boundary_conditions" in vessel:
            if location in vessel["boundary_conditions"]:
                ids_of_cap_vessels.append(vessel["vessel_id"])
    return ids_of_cap_vessels

