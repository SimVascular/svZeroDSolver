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
from .blocks import wire


def check_block_pair_flow_consistency(bA, bB):
    if bB.name not in bA.connecting_block_list:
        raise Exception("Block " + bB.name + " not in connecting list of " + bA.name)
    else:
        id_bB = bA.connecting_block_list.index(bB.name)

    if bA.name not in bB.connecting_block_list:
        raise Exception("Block " + bA.name + " not in connecting list of " + bB.name)
    else:
        id_bA = bB.connecting_block_list.index(bA.name)

    if bA.flow_directions[id_bB] * bB.flow_directions[id_bA] != -1:
        print("Flow direction of " + bB.name + " :", bB.flow_directions[id_bA])
        print("Flow direction of " + bA.name + " :", bA.flow_directions[id_bB])
        raise Exception(
            "Flow directions of " + bA.name + " donot conform to that of " + bB.name
        )


def connect_blocks_by_inblock_list(block_list):

    connectivity = []

    wire_dict = {}

    bnames = [_.name for _ in block_list]

    # Check if connection definition is consistent
    for bA in block_list:
        for bBnm in bA.connecting_block_list:
            bB = block_list[bnames.index(bBnm)]
            check_block_pair_flow_consistency(bA, bB)

    # If you reached here, it means each block has a consistent (connecting_block_list) and (flow_directions)
    for bA in block_list:
        i = -1
        id_bA = block_list.index(bA)
        for bBnm in bA.connecting_block_list:
            id_bB = bnames.index(bBnm)
            bB = block_list[id_bB]
            i += 1  # i is the index at which block, bB, is located in block bA's connecting_block_list
            if bA.flow_directions[i] == +1 and (id_bA, id_bB) not in connectivity:
                name_wire = bA.name + "_" + bB.name
                connecting_elements = (block_list[id_bA], block_list[id_bB])
                # wire_dict[name_wire] = wire(connecting_elements,name=name_wire)
                connectivity.append(
                    (id_bA, id_bB)
                )  # connectivity stores pair-wise tuples of indices of the blocks that are connected; basically, if block 1 is connected to block 2 and the flow goes from block 1 to block 2, then connectivity will store a 2-element tuple, where the first element is the index at which block 1 is stored in block_list and the 2nd element is the index at which block 2 is stored in block_list. if the flow goes from block 2 to block 1, then connectivity will store a 2-element tuple, where the first element is the index at which block 2 is stored in block_list and the 2nd element is the index at which block 1 is stored in block_list.
            elif bA.flow_directions[i] == -1:
                name_wire = bB.name + "_" + bA.name
                connecting_elements = (block_list[id_bB], block_list[id_bA])
            #     block_list[id_bA].add_connecting_wire(name_wire)
            #     block_list[id_bB].add_connecting_wire(name_wire)
            else:
                continue  # if this line is executed, then the next two lines (wire_dict[name_wire] = ... and block_list[id_bA] = ...) will not be executed
            wire_dict[name_wire] = wire(connecting_elements, name=name_wire)
            block_list[id_bA].add_connecting_wire(name_wire)

    return wire_dict


def check_block_connection(block):
    if len(block.flow_directions) != block.num_connections:
        print("Block name: " + block.name)
        print("Block number of flows: ", len(block.flow_directions))
        print("Block number of eqs: ", block.num_connections)

        raise Exception(
            "Number of connections donot match the number of inflows+outflows for this block"
        )

    # print block.connecting_wires_list
    reorder_inblock_connectivity(block)


# Reorder blocks to have connecting_block_list and connecting_wires_list arranged in ascending flow_directions
# This will give robustness to initial ordering during setup


def reorder_inblock_connectivity(block):
    indx = sorted(
        range(len(block.flow_directions)), key=lambda k: block.flow_directions[k]
    )

    block.flow_directions = [block.flow_directions[_] for _ in indx]
    block.connecting_wires_list = [block.connecting_wires_list[_] for _ in indx]
    block.connecting_block_list = [block.connecting_block_list[_] for _ in indx]


# Function to compute number of equations from blocks and wires
def compute_neq(block_list, wire_dict):
    neq = 0
    block_vars = 0
    for b in block_list:
        neq += b.neq
        block_vars += b.num_block_vars

    # print("Number of equations : ",neq)

    print(
        "Number of unknowns = ", 2 * len(wire_dict.values()) + block_vars
    )  # wire_dict.values() gives me an iterable or whatever whose length is the number of wires in wire_dict (number of wires in our model). then we multiply by 2, because each wire has 2 solution variables (P and Q).
    print(
        "Number of equations = ", neq
    )  # number of unknowns (solutionv variables) = 2*len(wire_dict.values()) + block_vars
    if 2 * len(wire_dict.values()) + block_vars != neq:
        print("Expected number of variables : ", 2 * len(wire_dict) + block_vars)
        print("Number of equations = ", neq)
        raise Exception("Mismatch between number of variables and equations")

    return neq


def initialize_solution_structures(neq):
    # Return y,ydot
    return np.zeros(neq), np.zeros(
        neq
    )  # recall that neq = number of solution variables = num of unknowns. thus, the global solution vector, y, should be of length neq


def assign_global_ids(
    block_list, wire_dict
):  # this function is where aekaansh assigns the global ids for the solution variables for the wires and blocks

    # Ordering of solution variables :
    # P0,Q0,P1,Q1,...,Pn,Qn, V1,V2,..,Vm # note that "V" stands for internal solution variable (of a block)
    # so the ordering of solution variables in the global vector of solution variables is: wire solutions first and then blocks' internal solutions

    i = 0  # i = the index at which a solution variable/unknown is stored in the global vector of solution variables/unknowns

    var_name_list = []

    # note that a solution ID = the index at which a solution variable is located in the global vector of solution variables

    for (
        w
    ) in (
        wire_dict.values()
    ):  # assign the wire solutions here (i.e. each wire has a P and Q solution. recall that each block, ie resistance block, has 2 associated wires and thus each block has 4 associated solutions (Pin, Qin, Pout, Qin). so here, we are assigning those solution ids in the global solution vector for those P and Q solutions
        # note that because wire_dict is a dictionary, it is unordered and basically, everytime we call wire_dict and loop through its values or keys or whatever, there is no set order of wires that we will follow and loop through.
        w.LPN_solution_ids = [i, i + 1]
        var_name_list.append("P_" + w.name)
        var_name_list.append("Q_" + w.name)
        i += 2

    for (
        b
    ) in (
        block_list
    ):  # here, we assign the solution ids for the internal solutions of the LPNBlocks
        b.LPN_solution_ids = []
        for j in range(b.num_block_vars):
            b.LPN_solution_ids.append(i)
            var_name_list.append("var_" + str(j) + "_" + b.name)
            i += 1

    offset = 0
    for b in block_list:
        for local_id in range(
            b.num_block_vars + 2 * len(b.connecting_block_list)
        ):  # note that b.num_block_vars+2*len(b.connecting_block_list) = the total number of solution variables/unknowns associated with this LPNBlock. len(b.connecting_block_list) is the number of wires (and blocks) attached to the current LPNBlock and this number is multiplied by 2 because each wire has 2 solutions (P and Q). then, the block also has internal solutions, where the number of internal solutions that it has is = b.num_block_vars
            b.global_col_id.append(
                b.eqids(wire_dict, local_id)
            )  # b.eqids returns the index at which the block's solution variable corresponding to local_id is located in the global vector of solution variables/unknowns.
        for local_id in range(b.neq):
            b.global_row_id += [offset + local_id]
        b.global_col_id = np.array(b.global_col_id)
        b.global_row_id = np.array(b.global_row_id)
        offset += b.neq
        meshgrid = np.array(np.meshgrid(b.global_row_id, b.global_col_id)).T.reshape(
            -1, 2
        )
        b.flat_row_ids, b.flat_col_ids = meshgrid[:, 0], meshgrid[:, 1]
        # recall that global_col_id is a list of the indices at which this LPNBlock's associated solution variables (Pin, Qin, Pout, Qout, and internal solutions) are stored in the global vector of solution variables/unknowns

    # print var_name_list

    return var_name_list
