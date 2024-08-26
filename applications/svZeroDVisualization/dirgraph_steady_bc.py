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

