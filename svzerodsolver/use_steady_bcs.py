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


class altered_bc_block:  # data structure to help restore the internal variables in the RCR and coronary BCs
    def __init__(self, block_name, bc_type, vessel_id, location, parameter_list):
        self.block_name = block_name
        self.bc_type = bc_type
        self.vessel_id = vessel_id
        self.location = location
        self.parameter_list = parameter_list  # parameters needed to compute the steady-state value of the internal variables


def compute_time_averaged_bc_value_for_single_cardiac_cycle(
    time, bc_values, cardiac_cycle_period
):
    """
    Inputs:
        np.array time
            = defined over a single cardiac cycle (time[-1] = time[0] + cardiac_cycle_period)
        np.array bc_values
            = defined over a single cardiac cycle (bc_values[0] == bc_values[-1])
        float cardiac_cycle_period
    """
    time_averaged_value = (1.0 / cardiac_cycle_period) * np.trapz(bc_values, time)
    return time_averaged_value


def get_bc_name_to_index_map(parameters):
    bc_name_to_index_map = {}
    for i in range(len(parameters["boundary_conditions"])):
        bc_name = parameters["boundary_conditions"][i]["bc_name"]
        bc_name_to_index_map[bc_name] = i
    return bc_name_to_index_map


def create_vessel_id_to_boundary_condition_map(parameters):
    vessel_id_to_boundary_condition_map = {}
    for vessel in parameters["vessels"]:
        if "boundary_conditions" in vessel:
            vessel_id = vessel["vessel_id"]
            vessel_id_to_boundary_condition_map[vessel_id] = {}
            for location, bc_name in vessel["boundary_conditions"].items():
                for boundary_condition in parameters["boundary_conditions"]:
                    if boundary_condition["bc_name"] == bc_name:
                        vessel_id_to_boundary_condition_map[vessel_id][
                            location
                        ] = boundary_condition
    parameters[
        "vessel_id_to_boundary_condition_map"
    ] = vessel_id_to_boundary_condition_map


def get_ids_of_cap_vessels(parameters, location):  # location == "inlet" or "outlet"
    ids_of_cap_vessels = []
    for vessel in parameters["vessels"]:
        if "boundary_conditions" in vessel:
            if location in vessel["boundary_conditions"]:
                ids_of_cap_vessels.append(vessel["vessel_id"])
    return ids_of_cap_vessels


def convert_unsteady_bcs_to_steady(parameters):
    """
    Convert unsteady BCs into equivalent steady BC by
        1. using the mean values of the pulsatile BCs as the values for the equivalent steady BCs
        2. removing capacitors from BCs (convert RCR and CORONARY into equivalent resistance-only BCs)

    Caveats:
        - removing the capacitors from the RCR and coronary BCs changes the 0d simulation results outputted by svzerodsolver.solver.run_network_util() because the internal variables originally present in the RCR and Coronary BCs are removed. These internal variables must be restored (and set to their steady state values) to use the steady state solutions as the initial conditions for the pulsatile simulations. The internal variables will be restored using the restore_internal_variables_for_capacitance_based_bcs() function.
    """
    # get helper data structures
    create_vessel_id_to_boundary_condition_map(parameters)
    vessel_id_to_boundary_condition_map = parameters[
        "vessel_id_to_boundary_condition_map"
    ]
    bc_name_to_index_map = get_bc_name_to_index_map(parameters)
    bc_identifiers = {"FLOW": "Q", "PRESSURE": "P", "CORONARY": "Pim"}

    # convert unsteady BCs to steady
    locations = ["inlet", "outlet"]
    altered_bc_blocks = (
        []
    )  # a list of the modified BC blocks (modified meaning "converted from RCR or coronary to resistance"); these BC blocks lost their internal variable in the conversion process
    for location in locations:
        ids_of_cap_vessels = get_ids_of_cap_vessels(parameters, location)
        for vessel_id in ids_of_cap_vessels:
            bc_name = vessel_id_to_boundary_condition_map[vessel_id][location][
                "bc_name"
            ]
            bc_type = vessel_id_to_boundary_condition_map[vessel_id][location][
                "bc_type"
            ]
            index_of_bc = bc_name_to_index_map[bc_name]
            if bc_type in bc_identifiers:
                # convert pulsatile BCs into steady BCs, whose value is the mean value of the pulsatile waveform
                bc_identifier = bc_identifiers[bc_type]
                time = vessel_id_to_boundary_condition_map[vessel_id][location][
                    "bc_values"
                ]["t"]
                bc_values = vessel_id_to_boundary_condition_map[vessel_id][location][
                    "bc_values"
                ][bc_identifier]
                if len(time) < 2:
                    message = (
                        "Error. Boundary condition, "
                        + bc_name
                        + ", must be prescribed for at least 2 time points."
                    )
                    raise RuntimeError(message)
                cardiac_cycle_period = time[-1] - time[0]
                time_averaged_value = (
                    compute_time_averaged_bc_value_for_single_cardiac_cycle(
                        time, bc_values, cardiac_cycle_period
                    )
                )
                parameters["boundary_conditions"][index_of_bc]["bc_values"]["t"] = [
                    time[0],
                    time[-1],
                ]
                parameters["boundary_conditions"][index_of_bc]["bc_values"][
                    bc_identifier
                ] = [time_averaged_value, time_averaged_value]
            elif bc_type == "RCR":
                # remove capacitor from RCR BC
                Rp = parameters["boundary_conditions"][index_of_bc]["bc_values"]["Rp"]
                Rd = parameters["boundary_conditions"][index_of_bc]["bc_values"]["Rd"]
                Pd = parameters["boundary_conditions"][index_of_bc]["bc_values"]["Pd"]
                equivalent_R = Rp + Rd  # add resistances in series
                parameters["boundary_conditions"][index_of_bc]["bc_values"] = {
                    "R": equivalent_R,
                    "Pd": Pd,
                }
                parameters["boundary_conditions"][index_of_bc]["bc_type"] = "RESISTANCE"
                altered_bc_blocks.append(
                    altered_bc_block(
                        "BC" + str(vessel_id) + "_" + location,
                        bc_type,
                        vessel_id,
                        location,
                        [Rp],
                    )
                )
            if bc_type == "CORONARY":
                # remove capacitors from coronary BC
                Ra1 = parameters["boundary_conditions"][index_of_bc]["bc_values"]["Ra1"]
                Ra2 = parameters["boundary_conditions"][index_of_bc]["bc_values"]["Ra2"]
                Rv1 = parameters["boundary_conditions"][index_of_bc]["bc_values"]["Rv1"]
                Cc = parameters["boundary_conditions"][index_of_bc]["bc_values"]["Cc"]
                Pim = parameters["boundary_conditions"][index_of_bc]["bc_values"]["Pim"]
                Pv_distal_pressure = parameters["boundary_conditions"][index_of_bc][
                    "bc_values"
                ]["P_v"]
                equivalent_R = Ra1 + Ra2 + Rv1  # add resistances in series
                parameters["boundary_conditions"][index_of_bc]["bc_values"] = {
                    "R": equivalent_R,
                    "Pd": Pv_distal_pressure,
                }
                parameters["boundary_conditions"][index_of_bc]["bc_type"] = "RESISTANCE"
                altered_bc_blocks.append(
                    altered_bc_block(
                        "BC" + str(vessel_id) + "_" + location,
                        bc_type,
                        vessel_id,
                        location,
                        [Ra1, Ra2, Cc, Pim],
                    )
                )
    return parameters, altered_bc_blocks


def restore_internal_variables_for_capacitance_based_bcs(
    y_f, ydot_f, var_name_list_f, altered_bc_blocks
):
    """
    Restore the internal variables (and set them to their steady state values) for the boundary condition blocks with capacitors (RCR and coronary).
    """
    y0 = copy.deepcopy(y_f)
    ydot0 = copy.deepcopy(ydot_f)
    var_name_list = copy.deepcopy(var_name_list_f)
    for block in altered_bc_blocks:
        if block.bc_type in ["RCR", "CORONARY"]:
            if block.location == "outlet":
                if block.bc_type == "RCR":
                    # get BC parameters
                    Rp = block.parameter_list[0]

                    # get soltns
                    Pin = y0[
                        var_name_list.index(
                            "P_V" + str(block.vessel_id) + "_" + block.block_name
                        )
                    ]
                    Qin = y0[
                        var_name_list.index(
                            "Q_V" + str(block.vessel_id) + "_" + block.block_name
                        )
                    ]

                    # compute value of internal variable
                    P_internal = Pin - Rp * Qin
                    var_name_list.append("var_0_" + block.block_name)
                    y0 = np.append(y0, np.array(P_internal))
                elif block.bc_type == "CORONARY":
                    # get BC parameters
                    Ra1 = block.parameter_list[0]
                    Ra2 = block.parameter_list[1]
                    Cc = block.parameter_list[2]
                    Pim = block.parameter_list[3]

                    # get soltns
                    Pin = y0[
                        var_name_list.index(
                            "P_V" + str(block.vessel_id) + "_" + block.block_name
                        )
                    ]
                    Qin = y0[
                        var_name_list.index(
                            "Q_V" + str(block.vessel_id) + "_" + block.block_name
                        )
                    ]

                    # compute value of internal variable
                    Pd = Pin - Qin * (Ra1 + Ra2)
                    volume_internal = Cc * (Pd - Pim)
                    var_name_list.append("var_0_" + block.block_name)
                    y0 = np.append(y0, np.array(volume_internal))
                ydot0 = np.append(
                    ydot0, np.zeros(1)
                )  # the time derivative of the soltn is zero b/c at steady-state
            else:
                message = "Error. This function does not work for inlet RCR or inlet coronary BCs."
                raise RuntimeError(message)
    return y0, ydot0, var_name_list
