import re
import numpy as np

from collections import defaultdict


def get_vessel_id_to_length_map(parameters):
    vessel_id_to_length_map = {}
    for vessel in parameters["vessels"]:
        vessel_id = vessel["vessel_id"]
        vessel_length = vessel["vessel_length"]
        vessel_id_to_length_map[vessel_id] = vessel_length
    return vessel_id_to_length_map


def initialize_0d_results_dict_branch(parameters, zero_d_time):
    zero_d_results = {"flow": {}, "pressure": {}, "distance": {}}
    num_time_pts = len(zero_d_time)
    branch_segment_ids = defaultdict(list)  # {branch_id : [branch_segment_ids]}
    branch_id_to_vessel_id_map = defaultdict(
        dict
    )  # {branch_id : {branch_segment_id : vessel_id}}
    vessel_id_to_length_map = get_vessel_id_to_length_map(parameters)

    for vessel in parameters["vessels"]:
        vessel_id = vessel["vessel_id"]
        vessel_name = vessel["vessel_name"]
        vessel_name_split = vessel_name.split("_")
        branch_id = int(
            (re.match(r"([a-z]+)([0-9]+)", vessel_name_split[0], re.I)).groups()[1]
        )
        branch_segment_id = int(
            (re.match(r"([a-z]+)([0-9]+)", vessel_name_split[1], re.I)).groups()[1]
        )
        branch_segment_ids[branch_id].append(branch_segment_id)
        branch_id_to_vessel_id_map[branch_id][branch_segment_id] = vessel_id

    for branch_id in branch_segment_ids:
        num_nodes_for_branch = max(branch_segment_ids[branch_id]) + 2
        for qoi in zero_d_results:
            zero_d_results[qoi][branch_id] = np.zeros(
                (num_nodes_for_branch, num_time_pts)
            )
        zero_d_results["distance"][branch_id] = np.zeros(num_nodes_for_branch)
        branch_segment_ids[
            branch_id
        ].sort()  # sort by ascending order of branch_segment_id
        counter = 0
        for branch_segment_id in branch_segment_ids[branch_id]:
            vessel_id = branch_id_to_vessel_id_map[branch_id][branch_segment_id]
            zero_d_results["distance"][branch_id][branch_segment_id + 1] = (
                zero_d_results["distance"][branch_id][counter]
                + vessel_id_to_length_map[vessel_id]
            )
            counter += 1
    zero_d_results["time"] = zero_d_time
    return zero_d_results


def get_vessel_id_to_vessel_name_map(parameters):
    vessel_id_to_vessel_name_map = {}
    for vessel in parameters["vessels"]:
        vessel_id = vessel["vessel_id"]
        vessel_name = vessel["vessel_name"]
        vessel_id_to_vessel_name_map[vessel_id] = vessel_name
    return vessel_id_to_vessel_name_map


def reformat_network_util_results_branch(
    zero_d_time, results_0d, var_name_list, parameters
):
    """
    Purpose:
        Reformat the 0d simulation results for just the branches into a dictionary (zero_d_results)
    Inputs:
        np.array zero_d_time
            = np.array of simulated time points
        np.array results_0d
            = np.array of the 0d simulation results, where the rows correspond to the each simulated time point and column j corresponds to the 0d solution for the solution variable name in var_name_list[j]
        list var_name_list
            = list of the 0d simulation results' solution variable names; most of the items in var_name_list are the QoIs + the names of the wires used in the 0d model (the wires connecting the 0d LPNBlock objects), where the wire names are usually comprised of the wire's inlet block name + "_" + the wire's outlet block name
                Example:
                    for var_name_list = ['P_V6_BC6_outlet', 'Q_V6_BC6_outlet'], then results_0d[:, i] holds the pressure (i = 0) or flow rate simulation result (i = 1) (both as np.arrays) for wire R6_BC6_outlet. This wire connects a resistance vessel block to an outlet BC block (specifically for vessel segment #6)
        dict parameters
            -- created from function utils.extract_info_from_solver_input_file
    Returns:
        dict zero_d_results
            =   {
                    "time" : 1d np.array of simulated time points,

                    "distance" : {branch_id : 1d np.array of distance of the branch's 0d nodes along the centerline}

                    "flow" : {branch_id : 2d np.array of flow rate where each row represents a 0d node on the branch and each column represents a time point,

                    "pressure" : {branch_id : 2d np.array of pressure where each row represents a 0d node on the branch and each column represents a time point}
                }

            - examples:
                1. plt.plot(zero_d_results["time"], zero_d_results["flow"][branch_id][0, :])
                        --> plot time vs the 0d flow waveform for the 0th node on branch, branch_id; this yields a plot that shows how the flow rate changes over time

                2. plt.plot(zero_d_results["distance"], zero_d_results["pressure"][branch_id][:, -1])
                        --> plot centerline distance vs the 0d pressure (at the last simulated time step); this yields a plot that shows how the pressure changes along the axial dimension of a vessel
    """
    qoi_map = {"Q": "flow", "P": "pressure"}
    zero_d_results = initialize_0d_results_dict_branch(parameters, zero_d_time)
    vessel_id_to_vessel_name_map = get_vessel_id_to_vessel_name_map(parameters)

    for i in range(len(var_name_list)):
        if "var" not in var_name_list[i]:  # var_name_list[i] == wire_name
            # the only possible combination of wire connections are: 1) vessel <--> vessel 2) vessel <--> junction 3) vessel <--> boundary condition; in all of these cases, there is at least one vessel block ("V") in each wire_name
            if "V" not in var_name_list[i]:
                message = "Error. It is expected that every wire in the 0d model must be connected to at least one vessel block."
                raise RuntimeError(message)
            else:
                # parse var_name
                var_name_split = var_name_list[i].split("_")
                qoi_header = var_name_split[0]

                if var_name_split[1].startswith(
                    "V"
                ):  # the wire connected downstream of this vessel block
                    vessel_id = int(var_name_split[1][1:])
                    vessel_name = vessel_id_to_vessel_name_map[vessel_id]
                    vessel_name_split = vessel_name.split("_")
                    branch_id = int(
                        (
                            re.match(r"([a-z]+)([0-9]+)", vessel_name_split[0], re.I)
                        ).groups()[1]
                    )
                    branch_segment_id = int(
                        (
                            re.match(r"([a-z]+)([0-9]+)", vessel_name_split[1], re.I)
                        ).groups()[1]
                    )
                    branch_node_id = branch_segment_id + 1
                    zero_d_results[qoi_map[qoi_header]][branch_id][
                        branch_node_id, :
                    ] = results_0d[:, i]
                else:  # need to find the inlet wire/node of the branch
                    if var_name_split[1].startswith(
                        "BC"
                    ):  # inlet wire/node of the branch
                        vessel_id = int(var_name_split[3][1:])
                        vessel_name = vessel_id_to_vessel_name_map[vessel_id]
                        vessel_name_split = vessel_name.split("_")
                        branch_id = int(
                            (
                                re.match(
                                    r"([a-z]+)([0-9]+)", vessel_name_split[0], re.I
                                )
                            ).groups()[1]
                        )
                        branch_segment_id = int(
                            (
                                re.match(
                                    r"([a-z]+)([0-9]+)", vessel_name_split[1], re.I
                                )
                            ).groups()[1]
                        )
                        if branch_segment_id != 0:
                            message = "Error. branch_segment_id should be 0 here because we are at the inlet wire of the branch."
                            raise RuntimeError(message)
                        else:
                            branch_node_id = branch_segment_id
                            zero_d_results[qoi_map[qoi_header]][branch_id][
                                branch_node_id, :
                            ] = results_0d[:, i]
                    elif var_name_split[1].startswith(
                        "J"
                    ):  # this wire could either be 1) the inlet wire/node of a branch or 2) some internal wire in the branch (where that internal wire is 1) connecting 2 vessel blocks or 2) connecting a vessel block and a junction block) (and we dont care about internal wires)
                        vessel_id = int(var_name_split[2][1:])
                        vessel_name = vessel_id_to_vessel_name_map[vessel_id]
                        vessel_name_split = vessel_name.split("_")
                        branch_id = int(
                            (
                                re.match(
                                    r"([a-z]+)([0-9]+)", vessel_name_split[0], re.I
                                )
                            ).groups()[1]
                        )
                        branch_segment_id = int(
                            (
                                re.match(
                                    r"([a-z]+)([0-9]+)", vessel_name_split[1], re.I
                                )
                            ).groups()[1]
                        )
                        if (
                            branch_segment_id == 0
                        ):  # this is the inlet wire/node of the branch
                            branch_node_id = branch_segment_id
                            zero_d_results[qoi_map[qoi_header]][branch_id][
                                branch_node_id, :
                            ] = results_0d[:, i]
                        else:  # this is an internal wire/node in the branch
                            pass  # do nothing here, since we are ignoring the internal wires where the vessel block is connected downstream of the wire (and the junction block is connected upstream of the wire)
                    else:
                        message = 'Error. It is not possible for a block name to begin with something other than, "V", "J", or "BC".'
                        raise RuntimeError(message)

    return zero_d_results
