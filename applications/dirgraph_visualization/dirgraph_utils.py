#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code will create a directed graph by reading a .json file.
"""
import dirgraph_wire as ntwku
import dirgraph_steady_bc
import dirgraph_connections
import numpy as np
import pandas as pd
import json
import pdb
import sys
import re
import os
import copy

try:
    from tqdm import tqdm
except ImportError:
    pass

try:
    import matplotlib.pyplot as plt  # only needed if you want to visualize the 0d model as a directed graph
except ImportError:
    print(
        "\nmatplotlib.pyplot not found. matplotlib.pyplot is needed only if you want to visualize your 0d model as a directed graph.")

try:
    import networkx as nx  # only needed if you want to visualize the 0d model as a directed graph
except ImportError:
    print("\nnetworkx not found. networkx is needed only if you want to visualize your 0d model as a directed graph.")

try:
    from profilehooks import profile  # only needed if you want to profile this script
except ImportError:
    print("\nprofilehooks.profile not found. profilehooks.profile is needed only if you want to profile this script.")

try:
    import json
except ImportError:
    print("\njson not found.")

import importlib
import argparse
from collections import defaultdict

# Loads a json impact file and extracts necessary information to draw an electrical circuit
def load_json_input_file(fpath, name_type, inlet_block_type):
    with open(fpath, 'rb') as fp:
        d = json.load(fp)
    blocks = {}  # {block_name : block_object}
    d.update({"blocks": blocks})
    dirgraph_steady_bc.create_block_to_boundary_condition_map(d)
    d.update({"block_names": list(d["blocks"].keys())})

    if 'vessels' not in d or len(d['vessels']) == 0:
        df_vessels = pd.DataFrame()

    else:
        df_vessels = pd.DataFrame(dict(
        inlet=[x.get('boundary_conditions', {}).get('inlet') for x in d['vessels']],
        outlet=[x.get('boundary_conditions', {}).get('outlet') for x in d['vessels']],
        name=[x['vessel_name'] for x in d['vessels']],
        vessel_id=[x['vessel_id'] for x in d['vessels']]
    ))


    if 'junctions' not in d:
        df_junctions = pd.DataFrame()
        df_junctions_expanded = pd.DataFrame()

    else:
        if inlet_block_type:
            inlet = 'inlet_blocks'
            outlet = 'outlet_blocks'
        else:
            inlet = 'inlet_vessels'
            outlet = 'outlet_vessels'

            # Create separate lists for inlets and outlets
        junction_inlets = []
        junction_outlets = []
        for junction in d['junctions']:
            junction_name = junction['junction_name']
            for block in junction.get(inlet, []):
                junction_inlets.append({'junction_name': junction_name, 'block_name': "V" + str(block), 'direction': 'inlet'})
            for block in junction.get(outlet, []):
                junction_outlets.append({'junction_name': junction_name, 'block_name': "V" + str(block), 'direction': 'outlet'})

        # Create DataFrames from the lists
        df_junctions_inlets = pd.DataFrame(junction_inlets)
        df_junctions_outlets = pd.DataFrame(junction_outlets)

        # Concatenate the DataFrames to form a unified junctions DataFrame
        df_junctions = pd.concat([df_junctions_inlets, df_junctions_outlets], ignore_index=True)

        df_junctions_expanded = pd.concat([df_junctions_inlets, df_junctions_outlets], ignore_index=True)

        create_junction_blocks(d, df_junctions_expanded, name_type)

    bcs = d['boundary_conditions']
    df_bcs = pd.DataFrame(dict(
        name=[x['bc_name'] for x in bcs],
        type=[x['bc_type'] for x in bcs],
    ))



    if 'valves' not in d:
        df_valves = pd.DataFrame(columns=['name', 'type', 'upstream', 'downstream'])
    else:
        valves = d['valves']
        df_valves = pd.DataFrame(dict(
            name = [x['name'] for x in valves],
            type = [x['type'] for x in valves],
            upstream = [x['params']['upstream_block'] for x in valves],
            downstream = [x['params']['downstream_block'] for x in valves]
        ))
        create_valve_blocks(d, df_vessels, df_valves, name_type)


    if 'chambers' not in d:
        df_chambers= pd.DataFrame(columns=['name', 'type'])

    else:
        chambers = d['chambers']
        df_chambers = pd.DataFrame(dict(
            name = [x['name'] for x in chambers],
            type=[x['type'] for x in chambers],
        ))

        create_chamber_blocks(d, df_chambers, df_junctions_expanded, df_valves)

    create_outlet_bc_blocks(d, df_valves, name_type)
    create_inlet_bc_blocks(d, df_valves, name_type)
    create_vessel_blocks(d, df_vessels, df_junctions_expanded, df_valves, name_type)
    return d



def create_valve_blocks(d, df_vessels, df_valves, name_type):
    valve_blocks = {}  # {block_name: block_object}

    if name_type == 'id':
        vessel_id_map = get_vessel_name_to_vessel_id_map(d)

    # vessel_list = get_vessel_list(d, 'name')
    bc_name_to_index_map = dirgraph_steady_bc.get_bc_name_to_index_map(d)

    def process_row(row):
        valve_name = row['name']
        ups_block = row['upstream']
        down_block = row['downstream']
        connecting_block_list = []
        flow_directions = []

        # Upstream Block Processing
        if ups_block in bc_name_to_index_map:
            bc_block_name = "BC" + valve_name + "_inlet"
            connecting_block_list.append(bc_block_name)
            flow_directions.append(-1)
        elif ups_block in df_vessels['name'].values:
            vessel_id_name = "V" + str(vessel_id_map[ups_block])
            connecting_block_list.append(vessel_id_name)
            flow_directions.append(-1)
        else:
            connecting_block_list.append(ups_block)
            flow_directions.append(-1)

        # Downstream Block Processing
        if down_block in bc_name_to_index_map:
            bc_block_name = "BC" + valve_name + "_outlet"
            connecting_block_list.append(bc_block_name)
            flow_directions.append(+1)
        elif down_block in df_vessels['name'].values:
            vessel_id_name = "V" + str(vessel_id_map[down_block])
            connecting_block_list.append(vessel_id_name)
            flow_directions.append(+1)
        else:
            connecting_block_list.append(down_block)
            flow_directions.append(+1)

        if (+1 in flow_directions) and (-1 in flow_directions):
            valve_blocks[valve_name] = ntwku.LPNBlock(
                connecting_block_list=connecting_block_list,
                name=valve_name,
                flow_directions=flow_directions
            )

    df_valves.apply(process_row, axis=1)
    d["blocks"].update(valve_blocks)



def create_chamber_blocks(d, df_chambers, df_junctions_expanded, df_valves):
    """
       Purpose:
           Create the chamber blocks for the 0d model.
       Inputs:
           dict parameters
               = created from function utils.extract_info_from_solver_input_file
           str name_type
               = str specified by either 'name' or 'id' that specifies whether
                 vessel names or ids should be used for each node
       Returns:
           void, but updates parameters["blocks"] to include the chamber, where
               parameters["blocks"] = {block_name : block_object}
       """
    chamber_blocks = {} # {block_name: block_object}

    def process_chamber(row):
        chamber_name = row['name']
        connecting_block_list = []
        flow_directions = []

        if not df_junctions_expanded.empty:
            # Junction Processing
            inlet_matches = df_junctions_expanded[(df_junctions_expanded['block_name'] == chamber_name) & (df_junctions_expanded['direction'] == 'inlet')]
            for _, match in inlet_matches.iterrows():
                connecting_block_list.append(match['junction_name'])
                flow_directions.append(+1)

            outlet_matches = df_junctions_expanded[(df_junctions_expanded['block_name'] == chamber_name) & (df_junctions_expanded['direction'] == 'outlet')]
            for _, match in outlet_matches.iterrows():
                connecting_block_list.append(match['junction_name'])
                flow_directions.append(-1)

        # Valve Processing
        valve_ups = df_valves[df_valves['upstream'] == chamber_name]
        for _, valve in valve_ups.iterrows():
            connecting_block_list.append(valve['name'])
            flow_directions.append(+1)

        valve_downs = df_valves[df_valves['downstream'] == chamber_name]
        for _, valve in valve_downs.iterrows():
            connecting_block_list.append(valve['name'])
            flow_directions.append(-1)

        if (+1 in flow_directions) and (-1 in flow_directions):
            chamber_blocks[chamber_name] = ntwku.LPNBlock(
                connecting_block_list=connecting_block_list,
                name=chamber_name,
                flow_directions=flow_directions
            )

    df_chambers.apply(process_chamber, axis=1)
    d["blocks"].update(chamber_blocks)


def create_junction_blocks(d, df_junctions_expanded, name_type):
    """
    Purpose:
        Create the junction blocks for the 0d model.
    Inputs:
        dict parameters
            = created from function utils.extract_info_from_solver_input_file
        str name_type
            = str specified by either 'name' or 'id' that specifies whether
              vessel names or ids should be used for each node
    Returns:
        void, but updates parameters["blocks"] to include the junction_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    if df_junctions_expanded.empty:
        return

    junction_blocks = {}  # {block_name: block_object}

    def process_junction(row):
        junction_name = row['junction_name']

        if not junction_name.startswith("J") or not junction_name[1:].isnumeric():
            message = f"Error. Joint name, {junction_name}, is not 'J' followed by numeric values. The 0D solver assumes that all joint names are 'J' followed by numeric values in the 0d solver input file. Note that the joint names are the same as the junction names."
            raise RuntimeError(message)

        connecting_block_list = []
        flow_directions = []

        # Inlet Processing
        inlet_matches = df_junctions_expanded[
            (df_junctions_expanded['junction_name'] == junction_name) & (df_junctions_expanded['direction'] == 'inlet')]
        for _, match in inlet_matches.iterrows():
            connecting_block_list.append(match['block_name'])
            flow_directions.append(-1)

        # Outlet Processing
        outlet_matches = df_junctions_expanded[(df_junctions_expanded['junction_name'] == junction_name) & (
                    df_junctions_expanded['direction'] == 'outlet')]
        for _, match in outlet_matches.iterrows():
            connecting_block_list.append(match['block_name'])
            flow_directions.append(+1)

        if (+1 in flow_directions) and (-1 in flow_directions):
            junction_blocks[junction_name] = ntwku.LPNBlock(
                connecting_block_list=connecting_block_list,
                name=junction_name,
                flow_directions=flow_directions
            )
        else:
            message = f"Error. Junction block, {junction_name}, must have at least 1 inlet connection and 1 outlet connection."
            raise RuntimeError(message)

    unique_junction_names = df_junctions_expanded['junction_name'].unique()
    for j_name in unique_junction_names:
        process_junction(df_junctions_expanded[df_junctions_expanded['junction_name'] == j_name].iloc[0])

    d["blocks"].update(junction_blocks)


def get_vessel_list(parameters, name_type):
    vessel_id_list = []
    for vessel in parameters["vessels"]:
        if name_type == 'id':
            vessel_id_list.append(vessel["vessel_id"])
        else:
            vessel_id_list.append(vessel['vessel_name'])
    return vessel_id_list


def get_vessel_block_helpers(d, df_vessels, df_junctions_expanded, df_valves, name_type):
    """
    Purpose:
        Create helper dictionaries to support the creation of the vessel blocks.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
        str name_type
            = str specified by either 'name' or 'id' that specifies whether
              vessel names or ids should be used for each node
    Returns:
        dict vessel_blocks_connecting_block_lists
            = {vessel_id : connecting_block_list}
        dict vessel_blocks_flow_directions
            = {vessel_id : flow_directions}
        dict vessel_blocks_names
            = {vessel_id : block_name}
    """
    name_type = name_type.lower()

    vessel_id_map = get_vessel_name_to_vessel_id_map(d)
    vessel_name_get = get_vessel_id_to_vessel_name_map(d) if name_type == 'name' else None

    # Initialize connecting block lists and flow directions
    df_vessels['connecting_block_list'] = [[] for _ in range(len(df_vessels))]
    df_vessels['flow_directions'] = [[] for _ in range(len(df_vessels))]
    df_vessels['block_name'] = df_vessels['vessel_id'].apply(
        lambda x: vessel_name_get[x] if name_type == 'name' else f"V{x}")

    # Process cap vessels
    for location in ["inlet", "outlet"]:
        ids_of_cap_vessels = dirgraph_steady_bc.get_ids_of_cap_vessels(d, location)
        direction = -1 if location == "inlet" else 1
        bc_block_names = [f"BC{v_id}_{location}" for v_id in ids_of_cap_vessels]

        cap_vessels_df = pd.DataFrame({
            'vessel_id': ids_of_cap_vessels,
            'bc_block_name': bc_block_names,
            'flow_direction': [direction] * len(ids_of_cap_vessels)
        })

        df_vessels = df_vessels.merge(cap_vessels_df, on='vessel_id', how='left')
        df_vessels['connecting_block_list'] = df_vessels.apply(lambda row: row['connecting_block_list'] + (
            [row['bc_block_name']] if pd.notna(row['bc_block_name']) else []), axis=1)
        df_vessels['flow_directions'] = df_vessels.apply(
            lambda row: row['flow_directions'] + ([row['flow_direction']] if pd.notna(row['flow_direction']) else []),
            axis=1)

        df_vessels.drop(columns=['bc_block_name', 'flow_direction'], inplace=True)


    for _, row in df_junctions_expanded.iterrows():
        vessel_id = row['block_name']
        id = int(vessel_id[1:])
        junction_name = row['junction_name']
        direction = +1 if row['direction'] == 'inlet' else -1

        df_vessels.loc[df_vessels['vessel_id'] == id, 'connecting_block_list'].apply(lambda x: x.append(junction_name))
        df_vessels.loc[df_vessels['vessel_id'] == id, 'flow_directions'].apply(lambda x: x.append(direction))

    # Process valves
    for _, row in df_valves.iterrows():
        blocks = [row['upstream'], row['downstream']]
        valve_name = row['name']
        for i, block in enumerate(blocks):
            if block in df_vessels['name'].values:
                vessel_id = df_vessels.loc[df_vessels['name'] == block, 'vessel_id'].values[0]
                df_vessels.loc[df_vessels['vessel_id'] == vessel_id, 'connecting_block_list'].apply(
                    lambda x: x.append(valve_name))
                flow_direction = +1 if i == 0 else -1
                df_vessels.loc[df_vessels['vessel_id'] == vessel_id, 'flow_directions'].apply(
                    lambda x: x.append(flow_direction))

    vessel_blocks_connecting_block_lists = df_vessels.set_index('vessel_id')['connecting_block_list'].to_dict()
    vessel_blocks_flow_directions = df_vessels.set_index('vessel_id')['flow_directions'].to_dict()
    vessel_blocks_names = df_vessels.set_index('vessel_id')['block_name'].to_dict()
    return vessel_blocks_connecting_block_lists, vessel_blocks_flow_directions, vessel_blocks_names



def create_vessel_blocks(d, df_vessels, df_junctions_expanded, df_valves, name_type):
    """
    Purpose:
        Create the vessel blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
        module custom_0d_elements_arguments
            = module to call custom 0d element arguments from
        str name_type
            = str specified by either 'name' or 'id' that specifies whether
              vessel names or ids should be used for each node
    Returns:
        void, but updates parameters["blocks"] to include the vessel_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    name_type = name_type.lower()
    vessel_blocks = {}  # {block_name : block_object}
    vessel_blocks_connecting_block_lists, vessel_blocks_flow_directions, vessel_blocks_names = get_vessel_block_helpers(
        d, df_vessels, df_junctions_expanded, df_valves, name_type)

    for vessel_id, block_name in vessel_blocks_names.items():
        connecting_block_list = vessel_blocks_connecting_block_lists[vessel_id]
        flow_directions = vessel_blocks_flow_directions[vessel_id]
        vessel_blocks[block_name] = ntwku.LPNBlock(
            connecting_block_list=connecting_block_list,
            name=block_name,
            flow_directions=flow_directions
        )

    d["blocks"].update(vessel_blocks)


def create_outlet_bc_blocks(d, df_valves, name_type):
    """
    Purpose:
        Create the outlet bc (boundary condition) blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
        module custom_0d_elements_arguments
            = module to call custom 0d element arguments from
        str name_type
            = str specified by either 'name' or 'id' that specifies whether
              vessel names or ids should be used for each node
    Returns:
        void, but updates parameters["blocks"] to include the outlet_bc_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    def process_boundary_condition(block_name, connecting_block_list, flow_directions):
        outlet_bc_blocks[block_name] = ntwku.LPNBlock(connecting_block_list=connecting_block_list, name=block_name,
                                                      flow_directions=flow_directions)

    name_type = name_type.lower()
    outlet_bc_blocks = {}
    outlet_vessels_of_model = dirgraph_steady_bc.get_ids_of_cap_vessels(d, "outlet")
    block_to_boundary_condition_map = d["block_to_boundary_condition_map"]
    vessel_name_get = get_vessel_id_to_vessel_name_map(d)

    for vessel_id in outlet_vessels_of_model:
        block_name = f"BC{vessel_id}_outlet"
        connecting_block_list = [vessel_name_get[vessel_id]] if name_type == 'name' else [f"V{vessel_id}"]
        flow_directions = [-1]
        process_boundary_condition(block_name,connecting_block_list, flow_directions)

    for _, valve in df_valves.iterrows():
        valve_name = valve['name']
        if valve_name not in block_to_boundary_condition_map or 'outlet' not in block_to_boundary_condition_map[
            valve_name]:
            continue
        block_name = f"BC{valve_name}_outlet"
        connecting_block_list = [valve_name]
        flow_directions = [-1]
        process_boundary_condition(block_name, connecting_block_list, flow_directions)

    d["blocks"].update(outlet_bc_blocks)

def create_inlet_bc_blocks(d, df_valves, name_type):
    """
    Purpose:
        Create the inlet bc (boundary condition) blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
        module custom_0d_elements_arguments
            = module to call custom 0d element arguments from
    Returns:
        void, but updates parameters["blocks"] to include the inlet_bc_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """

    name_type = name_type.lower()
    inlet_bc_blocks = {}
    inlet_vessels_of_model = dirgraph_steady_bc.get_ids_of_cap_vessels(d, "inlet")
    block_to_boundary_condition_map = d["block_to_boundary_condition_map"]
    vessel_name_get = get_vessel_id_to_vessel_name_map(d)

    # Process valves

    for _, valve in df_valves.iterrows():
        valve_name = valve['name']
        if valve_name not in block_to_boundary_condition_map or 'inlet' not in block_to_boundary_condition_map[
            valve_name]:
            continue
        block_name = f"BC{valve_name}_inlet"
        connecting_block_list = [valve_name]
        flow_directions = [+1]
        inlet_bc_blocks[block_name] = ntwku.LPNBlock(
            connecting_block_list=connecting_block_list,
            name=block_name,
            flow_directions=flow_directions
        )

    # Process vessels
    for vessel_id in inlet_vessels_of_model:
        block_name = f"BC{vessel_id}_inlet"
        connecting_block_list = [vessel_name_get[vessel_id]] if name_type == 'name' else [f"V{vessel_id}"]
        flow_directions = [+1]
        location = "inlet"
        inlet_bc_blocks[block_name] = ntwku.LPNBlock(
            connecting_block_list=connecting_block_list,
            name=block_name,
            flow_directions=flow_directions
        )

    d["blocks"].update(inlet_bc_blocks)


# @profile
def run_network_util(zero_d_solver_input_file_path, d, draw_directed_graph, output_dir):
    """
    Purpose:
        Run functions from network_util_NR to execute the 0d simulation and generate simulation results (pressure, flow rates).
    Inputs:
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
        boolean draw_directed_graph
            = True to visualize the 0d model as a directed graph using networkx -- saves the graph to a .png file (hierarchical graph layout) and a networkx .dot file; False, otherwise. .dot file can be opened with neato from graphviz to visualize the directed in a different format.
    Returns:
        list block_list
            = list of blocks (nodes)
        list connect_list
            = list of connections between blocks (edges)
    """
    block_list = list(d["blocks"].values())
    connect_list, wire_dict = dirgraph_connections.connect_blocks_by_inblock_list(block_list)
    case_name = zero_d_solver_input_file_path.rsplit('/', 1)[-1]
    # directed_graph_file_path = output_dir + '/' + os.path.splitext(case_name)[0] + "_directed_graph"
    directed_graph_file_path = os.path.join(output_dir, os.path.splitext(case_name)[0] + "_directed_graph")
    save_directed_graph(block_list, connect_list, directed_graph_file_path, draw_directed_graph)

    return block_list, connect_list


def get_vessel_id_to_length_map(parameters):
    vessel_id_to_length_map = {}
    for vessel in parameters["vessels"]:
        vessel_id = vessel["vessel_id"]
        vessel_length = vessel["vessel_length"]
        vessel_id_to_length_map[vessel_id] = vessel_length
    return vessel_id_to_length_map


def get_vessel_id_to_vessel_name_map(parameters):
    vessel_id_to_vessel_name_map = {}
    for vessel in parameters["vessels"]:
        vessel_id = vessel["vessel_id"]
        vessel_name = vessel["vessel_name"]
        vessel_id_to_vessel_name_map[vessel_id] = vessel_name
    return vessel_id_to_vessel_name_map

def get_vessel_name_to_vessel_id_map(parameters):
    vessel_name_to_vessel_id_map = {}
    for vessel in parameters["vessels"]:
        vessel_name = vessel["vessel_name"]
        vessel_id = vessel["vessel_id"]
        vessel_name_to_vessel_id_map[vessel_name] = vessel_id
    return vessel_name_to_vessel_id_map


def save_directed_graph(block_list, connect_list, directed_graph_file_path, draw_directed_graph):
    """
    Purpose:
        Visualize the 0d model as a directed graph -- save the graph in a hierarchical graph layout to a .png file; also save a networkx .dot file that can be opened with neato via graphviz to visualize the graph in a different layout.
    Inputs:
        list block_list
            = [list of all of the 0d LPNBlock objects]
        list connect_list
            = [list of (blockA_index, blockB_index)]
                where blockA and blockB are connected to each other, and blockA_index and blockB_index are the index locations at which the blockA and blockB objects are stored in block_list
        string directed_graph_file_path
            = name of the hierarchical graph .png file and networkx .dot file that will be saved
    Returns:
        void, but saves a .png file visualizing the 0d model as a directed graph, as well as a networkx .dot file that can be opened with neato via graphviz to visualize the graph in a different layout
    """
    G = nx.DiGraph()
    G.add_edges_from([(block_list[tpl[0]].name, block_list[tpl[1]].name) for tpl in connect_list])
    pos = nx.nx_pydot.pydot_layout(G,
                                   prog='dot')
    options = {"alpha": 0.2}
    nx.draw_networkx_nodes(G, pos, node_color='blue', node_size=300, **options)
    nx.draw_networkx_labels(G, pos, font_size='6')
    nx.draw_networkx_edges(G, pos, node_size=300)

    nx.nx_pydot.write_dot(G, directed_graph_file_path + ".dot")

    if draw_directed_graph == True:
        plt.figure(figsize=(20, 11))
        plt.tight_layout()
        plt.savefig(directed_graph_file_path + ".png", format="PNG", dpi=100)
        plt.show()
        plt.close("all")



def set_up_0d_network(zero_d_solver_input_file_path: object, output_dir, name_type: object, inlet_block: object = False, draw_directed_graph: object = False) -> object:
    """
    Purpose:
        Create all network_util_NR::LPNBlock objects for the 0d model and run the 0d simulation.
    Inputs:
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
        str name_type
            = str specified by either 'name' or 'id' that specifies whether
              vessel names or ids should be used for each node
        boolean draw_directed_graph
            = True to visualize the 0d model as a directed graph using networkx -- saves the graph to a .png file (hierarchical graph layout) and a networkx .dot file; False, otherwise. .dot file can be opened with neato from graphviz to visualize the directed in a different format.
        boolean use_custom_0d_elements
            = True to use user-defined, custom 0d elements in the 0d model; False, otherwire
        string custom_0d_elements_arguments_file_path
            = path to user-defined custom 0d element file
    Caveats:
        The save_results_branch option works only for 0d models with the branching structure where each vessel is modeled as a single branch with 1 or multiple sub-segments
    Returns:
        void
    """
    d = load_json_input_file(zero_d_solver_input_file_path, name_type, inlet_block)
    block_list, connect_list = run_network_util(
        zero_d_solver_input_file_path,
        d,
        draw_directed_graph=draw_directed_graph, output_dir=output_dir
    )
