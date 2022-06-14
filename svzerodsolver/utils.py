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

import os
from collections import OrderedDict
import json


def extract_info_from_solver_input_file(
    solver_input_file_path, one_d_inlet_segment_number=0, write_json=False
):
    """
    Purpose:
        Extract geometric and material properties prescribed in the 1d or 0d solver input file
    Inputs:
        string solver_input_file_path
            = name of the 1d or 0d solver input file
        int one_d_inlet_segment_number
            = inlet segment number of the 1d model (this number is usually zero)
    Returns:
        dict parameters = # see comments below for description of each item in parameters
            {
                "segment_names" : segment_names,
                "lengths" : lengths,
                "segment_0d_types" : segment_0d_types,
                "segment_0d_values" : segment_0d_values,
                "joint_name_to_segments_map" : joint_name_to_segments_map,
                "junction_types" : junction_types,
                "boundary_condition_types" : boundary_condition_types,
                "boundary_condition_datatable_names" : boundary_condition_datatable_names,
                "datatable_values" : datatable_values,
                "number_of_time_pts_per_cardiac_cycle" : number_of_time_pts_per_cardiac_cycle,
                "number_of_cardiac_cycles" : number_of_cardiac_cycles
            }
    """
    segment_0d_parameter_types = {
        "R": ["R"],
        "C": ["C"],
        "L": ["L"],
        "RC": ["R", "C"],
        "RL": ["R", "L"],
        "RCL": ["R", "C", "L"],
        "STENOSIS": ["R_poiseuille", "stenosis_coefficient"],
    }

    boundary_condition_parameter_types = {
        "FLOW": ["Q"],
        "PRESSURE": ["P"],
        "RESISTANCE": ["R", "Pd"],
        "RCR": ["Rp", "C", "Rd", "Pd"],
        "CORONARY": ["Ra1", "Ra2", "Ca", "Cc", "Rv1", "P_v", "Pim"],
    }

    segment_names = {}  # {segment number : segment name}
    lengths = {}  # {segment number : segment length}
    segment_0d_types = (
        {}
    )  # {segment number : 0d element type, i.e. R, RC, RCL, RL, L, C}
    segment_0d_values = (
        {}
    )  # {segment number : {"R" : val} or {"R" : val, "C" : val, "L : val"} }
    # example: if the 0d element type is R, then the 0d element value is [R_value]
    # example: if the 0d element type is RCL, then the 0d element values are [R_value, C_value, L_value]
    joint_node_numbers = {}  # {joint name : joint node number}
    joint_inlet_names = {}  # {joint name : joint inlet name}
    joint_outlet_names = {}  # {joint name : joint outlet name}
    joint_inlet_segments = (
        {}
    )  # {joint inlet name : [list of inlet segments at the joint]}
    joint_outlet_segments = (
        {}
    )  # {joint outlet name : [list of outlet segments at the joint]}
    joint_name_to_segments_map = (
        {}
    )  # { joint name : { "inlet" : [list of inlet segments at the joint], "outlet" : [list of outlet segments at the joint] } }
    junction_types = {}  # {joint name : junction type, i.e. NORMAL_JUNCTION}
    boundary_condition_types = {
        "inlet": {},
        "outlet": {},
    }  # {"inlet" : {segment number : type of boundary condition, i.e. "FLOW", "PRESSURE"}, "outlet" : {segment number : type of boundary condition, i.e. "RESISTANCE", "RCR", "CORONARY", "FLOW", "PRESSURE"}}
    boundary_condition_datatable_names = {
        "inlet": {},
        "outlet": {},
    }  # {"inlet" : {segment number : datatable name for boundary condition}, "outlet" : {segment number : corresponding datatable name for boundary condition}}
    datatable_names_to_boundary_condition_type_map = (
        {}
    )  # {datatable name for boundary condition : type of boundary condition, e.g. "RESISTANCE", "RCR"}
    datatable_values = (
        {}
    )  # {datatable name for boundary condition : [time1 value1 time2 value2 ...]}
    # value corresponds to the value of the BC, ie resistance value
    inlet_segments_of_model = (
        []
    )  # [model's inlet segments (these are attached to the inlet boundary conditions)]
    outlet_segments_of_model = (
        []
    )  # [model's outlet segments (these are attached to the outlet boundary conditions)]

    with open(solver_input_file_path) as f:
        while True:
            line = f.readline()
            if line == "":  # signifies end of file
                break
            else:
                if line.startswith("JOINT "):
                    line_list = line.split()
                    joint_name = line_list[1]
                    if (not joint_name.startswith("J")) or (
                        not joint_name[1].isnumeric()
                    ):
                        message = (
                            "Error in file, "
                            + solver_input_file_path
                            + ". Joint name, "
                            + joint_name
                            + ", is not 'J' followed by numeric values. The 0D solver assumes that all joint names are 'J' followed by numbers in the solver input file."
                        )
                        raise ValueError(message)
                    joint_node_numbers[joint_name] = int(line_list[2])
                    joint_inlet_names[joint_name] = line_list[3]
                    joint_outlet_names[joint_name] = line_list[4]
                elif line.startswith("JOINTINLET"):
                    line_list = line.split()
                    joint_in_name = line_list[1]
                    if (not joint_in_name.startswith("IN")) or (
                        not joint_in_name[2].isnumeric()
                    ):
                        message = (
                            "Error in file, "
                            + solver_input_file_path
                            + ". Joint inlet name, "
                            + joint_in_name
                            + ", is not 'IN' followed by numbers. This code assumes that all joint inlet names is 'IN' followed by numbers in the solver input file."
                        )
                        raise ValueError(message)
                    joint_inlet_segments[joint_in_name] = [
                        int(i) for i in line_list[3:]
                    ]  # inlet segment #'s at joint
                    if len(joint_inlet_segments[joint_in_name]) != int(line_list[2]):
                        message = (
                            "Error in file, "
                            + solver_input_file_path
                            + ". Number of inlet segments for joint inlet, "
                            + joint_in_name
                            + ", does not match the number of prescribed inlet segments."
                        )
                        raise ValueError(message)
                elif line.startswith("JOINTOUTLET"):
                    line_list = line.split()
                    joint_out_name = line_list[1]
                    if (not joint_out_name.startswith("OUT")) or (
                        not joint_out_name[3].isnumeric()
                    ):
                        message = (
                            "Error in file, "
                            + solver_input_file_path
                            + ". Joint outlet name, "
                            + joint_out_name
                            + ", is not 'OUT' followed by numbers. This code assumes that all joint outlet names is 'OUT' followed by numbers in the solver input file."
                        )
                        raise ValueError(message)
                    joint_outlet_segments[joint_out_name] = [
                        int(i) for i in line_list[3:]
                    ]  # outlet segment #'s at joint
                    if len(joint_outlet_segments[joint_out_name]) != int(line_list[2]):
                        message = (
                            "Error in file, "
                            + solver_input_file_path
                            + ". Number of outlet segments for joint outlet, "
                            + joint_out_name
                            + ", does not match the number of prescribed outlet segments."
                        )
                        raise ValueError(message)
                elif line.startswith(
                    "JUNCTION_MODEL"
                ):  # this section is only available in the 0d solver input file
                    line_list = line.split()
                    joint_name = line_list[1]
                    if (not joint_name.startswith("J")) or (
                        not joint_name[1].isnumeric()
                    ):
                        message = (
                            "Error in file, "
                            + solver_input_file_path
                            + ". Joint name, "
                            + joint_name
                            + ", is not 'J' followed by numeric values. The 0D solver assumes that all joint names are 'J' followed by numbers in the 0d solver input file. Note that the joint names are the same as the junction names."
                        )
                        raise ValueError(message)
                    junction_types[joint_name] = line_list[2]
                elif line.startswith(
                    "ELEMENT"
                ):  # this section is only available in the 0d solver input file
                    line_list = line.split()
                    segment_number = int(line_list[2])
                    segment_names[segment_number] = line_list[1]
                    lengths[segment_number] = float(line_list[4])
                    if line_list[13] != "NOBOUND":
                        boundary_condition_types["outlet"][segment_number] = line_list[
                            13
                        ]
                        boundary_condition_datatable_names["outlet"][
                            segment_number
                        ] = line_list[14]
                        outlet_segments_of_model.append(segment_number)

                        datatable_names_to_boundary_condition_type_map[
                            line_list[14]
                        ] = line_list[13]
                    segment_0d_types[segment_number] = line_list[15]
                    if line_list[15] in segment_0d_parameter_types:
                        segment_0d_values[segment_number] = {}
                        num_types = len(segment_0d_parameter_types[line_list[15]])
                        for i in range(num_types):
                            type = segment_0d_parameter_types[line_list[15]][i]
                            segment_0d_values[segment_number][type] = float(
                                line_list[16 + i]
                            )
                    else:
                        message = (
                            "Error in file, "
                            + solver_input_file_path
                            + ". 0D element type, "
                            + segment_0d_types[segment_number]
                            + ", is not recognized."
                        )
                        raise ValueError(message)
                elif line.startswith("DATATABLE"):
                    line_list = line.split()
                    datatable_name = line_list[1]
                    if line_list[2] == "LIST":
                        datatable_values[datatable_name] = OrderedDict()
                        bc_type = datatable_names_to_boundary_condition_type_map[
                            datatable_name
                        ]
                        if bc_type in boundary_condition_parameter_types:
                            num_param_types = len(
                                boundary_condition_parameter_types[bc_type]
                            )
                            for i in range(num_param_types):
                                line = f.readline()
                                line_list = line.split()
                                param_type = boundary_condition_parameter_types[
                                    bc_type
                                ][i]
                                if bc_type == "RCR" and param_type == "Pd":
                                    datatable_values[datatable_name][
                                        param_type
                                    ] = 0.0  # default distal pressure for RCR BC (since the DATATABLE for RCR BCs doesn't have a distal pressure listed); if the DATATABLE for the RCR BC did have a distal pressure listed (as the last value), then we would just use the command, datatable_values[datatable_name][param_type] = line_list[1], instead
                                else:
                                    datatable_values[datatable_name][
                                        param_type
                                    ] = float(line_list[1])
                                if bc_type in ["CORONARY", "FLOW", "PRESSURE"]:
                                    if i == num_param_types - 1:
                                        datatable_values[datatable_name]["t"] = []
                                        datatable_values[datatable_name][
                                            param_type
                                        ] = []

                                        while line_list[0] != "ENDDATATABLE":
                                            datatable_values[datatable_name][
                                                "t"
                                            ] = datatable_values[datatable_name][
                                                "t"
                                            ] + [
                                                float(line_list[0])
                                            ]
                                            datatable_values[datatable_name][
                                                param_type
                                            ] = datatable_values[datatable_name][
                                                param_type
                                            ] + [
                                                float(line_list[1])
                                            ]

                                            line = f.readline()
                                            line_list = line.split()
                        else:
                            message = (
                                "Boundary condition type, "
                                + bc_type
                                + ", is not recognized."
                            )
                            raise ValueError(message)
                    else:
                        message = (
                            "Error in file, "
                            + solver_input_file_path
                            + ". DATATABLE type, "
                            + line_list[2]
                            + ", is not allowed. Only allowed type is 'LIST'."
                        )
                        raise ValueError(message)
                elif line.startswith(
                    "INLET_BOUNDARY_CONDITION"
                ):  # this section is only available in the 0d solver input file
                    line_list = line.split()
                    segment_number = int(line_list[1])
                    inlet_segments_of_model.append(segment_number)
                    boundary_condition_types["inlet"][segment_number] = line_list[2]
                    boundary_condition_datatable_names["inlet"][
                        segment_number
                    ] = line_list[3]

                    datatable_names_to_boundary_condition_type_map[
                        line_list[3]
                    ] = line_list[2]

                elif line.startswith("SOLVEROPTIONS"):
                    line_list = line.split()
                    if line_list[0] == "SOLVEROPTIONS_0D":  # 0D solver options
                        number_of_time_pts_per_cardiac_cycle = int(line_list[1])
                        number_of_cardiac_cycles = int(line_list[2])
                    else:
                        message = (
                            "Error. Unidentified solver card, " + line_list[0] + "."
                        )
                        raise RuntimeError(message)

    segment_numbers_list = list(segment_names.keys())
    segment_numbers_list.sort()
    if len(segment_numbers_list) == 1:
        if segment_numbers_list != inlet_segments_of_model:
            message = (
                "Error in file, "
                + solver_input_file_path
                + ". This model has only 1 segment, segment #"
                + str(segment_numbers_list[0])
                + ", but this segment isn't equal to inlet_segments_of_model[0], "
                + str(inlet_segments_of_model[0])
                + "."
            )
            raise RuntimeError(message)
        elif segment_numbers_list != outlet_segments_of_model:
            message = (
                "Error in file, "
                + solver_input_file_path
                + ". This model has only 1 segment, segment #"
                + str(segment_numbers_list[0])
                + ", but this segment isn't equal to outlet_segments_of_model[0], "
                + str(outlet_segments_of_model[0])
                + "."
            )
            raise RuntimeError(message)

    for joint_name in list(joint_node_numbers.keys()):
        joint_name_to_segments_map[joint_name] = {
            "inlet": joint_inlet_segments[joint_inlet_names[joint_name]].copy(),
            "outlet": joint_outlet_segments[joint_outlet_names[joint_name]].copy(),
        }

    parameters = {"boundary_conditions": [], "vessels": [], "junctions": []}

    parameters["simulation_parameters"] = {
        "number_of_time_pts_per_cardiac_cycle": number_of_time_pts_per_cardiac_cycle,
        "number_of_cardiac_cycles": number_of_cardiac_cycles,
    }

    for bc_name, bc_type in datatable_names_to_boundary_condition_type_map.items():
        bc_values = datatable_values[bc_name]
        boundary_condition = {
            "bc_name": bc_name,
            "bc_type": bc_type,
            "bc_values": bc_values,
        }
        parameters["boundary_conditions"].append(boundary_condition)

    for vessel_id, vessel_name in segment_names.items():
        zero_d_element_type = segment_0d_types[vessel_id]
        zero_d_element_values = segment_0d_values[vessel_id]
        vessel_length = lengths[vessel_id]
        vessel = {
            "vessel_id": vessel_id,
            "vessel_name": vessel_name,
            "zero_d_element_type": zero_d_element_type,
            "zero_d_element_values": zero_d_element_values,
            "vessel_length": vessel_length,
        }

        boundary_conditions = None
        if (
            vessel_id in boundary_condition_datatable_names["inlet"]
            or vessel_id in boundary_condition_datatable_names["outlet"]
        ):
            boundary_conditions = {}
            for location in ["inlet", "outlet"]:
                if vessel_id in boundary_condition_datatable_names[location]:
                    boundary_conditions[location] = boundary_condition_datatable_names[
                        location
                    ][vessel_id]

        if boundary_conditions:
            vessel["boundary_conditions"] = boundary_conditions
        parameters["vessels"].append(vessel)

    for junction_name, junction_type in junction_types.items():
        inlet_vessels = joint_name_to_segments_map[junction_name]["inlet"]
        outlet_vessels = joint_name_to_segments_map[junction_name]["outlet"]
        junction = {
            "junction_name": junction_name,
            "junction_type": junction_type,
            "inlet_vessels": inlet_vessels,
            "outlet_vessels": outlet_vessels,
        }
        parameters["junctions"].append(junction)

    if write_json:
        json_file_name = os.path.splitext(solver_input_file_path)[0] + ".json"
        with open(json_file_name, "w") as outfile:
            json.dump(parameters, outfile, indent=4, sort_keys=True)
