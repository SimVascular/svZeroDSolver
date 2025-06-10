# coding=utf-8

# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
from dirgraph_wire import wire

def check_block_pair_flow_consistency(bA, bB):
    if bB.name not in bA.connecting_block_list:
        raise Exception('Block ' + bB.name + ' not in connecting list of ' + bA.name)
    else:
        id_bB = bA.connecting_block_list.index(bB.name)

    if bA.name not in bB.connecting_block_list:
        raise Exception('Block ' + bA.name + ' not in connecting list of ' + bB.name)
    else:
        id_bA = bB.connecting_block_list.index(bA.name)

    if bA.flow_directions[id_bB] * bB.flow_directions[id_bA] != -1:
        print('Flow direction of ' + bB.name + ' :', bB.flow_directions[id_bA])
        print('Flow direction of ' + bA.name + ' :', bA.flow_directions[id_bB])
        raise Exception('Flow directions of ' + bA.name + ' donot conform to that of ' + bB.name)


def connect_blocks_by_inblock_list(
        block_list):
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
                name_wire = bA.name + '_' + bB.name
                connecting_elements = (block_list[id_bA], block_list[id_bB])
                connectivity.append((id_bA,
                                     id_bB))  # connectivity stores pair-wise tuples of indices of the blocks that are connected; basically, if block 1 is connected to block 2 and the flow goes from block 1 to block 2, then connectivity will store a 2-element tuple, where the first element is the index at which block 1 is stored in block_list and the 2nd element is the index at which block 2 is stored in block_list. if the flow goes from block 2 to block 1, then connectivity will store a 2-element tuple, where the first element is the index at which block 2 is stored in block_list and the 2nd element is the index at which block 1 is stored in block_list.
            elif bA.flow_directions[i] == -1:
                name_wire = bB.name + '_' + bA.name
                connecting_elements = (block_list[id_bB], block_list[id_bA])
            else:
                continue  # if this line is executed, then the next two lines (wire_dict[name_wire] = ... and block_list[id_bA] = ...) will not be executed
            wire_dict[name_wire] = wire(connecting_elements, name=name_wire)
            block_list[id_bA].add_connecting_wire(name_wire)

    return connectivity, wire_dict


def connect_blocks_by_connectivity_list(block_list, connectivity):
    wire_dict = {}

    for e in connectivity:
        e1, e2 = e
        e1name = block_list[e1].name
        e2name = block_list[e2].name

        connecting_elements = (block_list[e1], block_list[e2])
        name_wire = e1name + '_' + e2name

        wire_dict[name_wire] = wire(connecting_elements, name=name_wire)

        if e2name not in block_list[e1].connecting_block_list:
            block_list[e1].add_connecting_wire(name_wire)
            block_list[e1].add_connecting_block(e2name, +1)

        if e1name not in block_list[e2].connecting_block_list:
            block_list[e2].add_connecting_wire(name_wire)
            block_list[e2].add_connecting_block(e1name, -1)
    return wire_dict


def check_block_connection(block):
    if len(block.flow_directions) != block.num_connections:
        print("Block name: " + block.name)
        print("Block number of flows: ", len(block.flow_directions))
        print("Block number of eqs: ", block.num_connections)

        raise Exception("Number of connections donot match the number of inflows+outflows for this block")

    # print block.connecting_wires_list
    reorder_inblock_connectivity(block)


# Reorder blocks to have connecting_block_list and connecting_wires_list arranged in ascending flow_directions
# This will give robustness to initial ordering during setup

def reorder_inblock_connectivity(block):
    indx = sorted(range(len(block.flow_directions)), key=lambda k: block.flow_directions[k])

    block.flow_directions = [block.flow_directions[_] for _ in indx]
    block.connecting_wires_list = [block.connecting_wires_list[_] for _ in indx]
    block.connecting_block_list = [block.connecting_block_list[_] for _ in indx]
