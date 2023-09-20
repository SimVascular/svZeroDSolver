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
"""This module holds the main execution routines of svZeroDSolver."""
import json
import pandas as pd
import numpy as np

import click

from .algebra import run_integrator
from .utils import (
    convert_unsteady_bcs_to_steady,
    create_blocks,
    format_results_to_dict,
    get_solver_params,
)


def run_from_config(parameters):
    """Run the svZeroDSolver.

    Args:
        config: Python dict of the configuration.

    Returns:
        Python dict with results.
    """

    y_initial = None
    ydot_initial = None
    sim_params = parameters["simulation_parameters"]
    if sim_params.get("steady_initial", True):
        steady_parameters = convert_unsteady_bcs_to_steady(parameters)

        # to run the 0d model with steady BCs to steady-state, simulate this model
        # with large time step size for an arbitrarily small number of cardiac cycles
        block_list, dofhandler = create_blocks(steady_parameters, steady=True)
        time_step_size, num_time_steps = get_solver_params(steady_parameters)
        (_, y_list, ydot_list) = run_integrator(
            block_list,
            dofhandler,
            num_time_steps,
            time_step_size,
            atol=sim_params.get("absolute_tolerance", 1e-8),
            max_iter=sim_params.get("maximum_nonlinear_iterations", 30),
        )
        y_initial = y_list[-1]
        ydot_initial = ydot_list[-1]

    block_list, dofhandler = create_blocks(parameters)
    time_step_size, num_time_steps = get_solver_params(parameters)
    (time_steps, y_list, ydot_list) = run_integrator(
        block_list,
        dofhandler,
        num_time_steps,
        time_step_size,
        y_initial,
        ydot_initial,
        atol=sim_params.get("absolute_tolerance", 1e-8),
        max_iter=sim_params.get("maximum_nonlinear_iterations", 30),
    )

    branch_result = format_results_to_dict(time_steps, y_list, block_list)

    result = pd.DataFrame(
        columns=[
            "name",
            "time",
            "flow_in",
            "flow_out",
            "pressure_in",
            "pressure_out",
        ]
    )
    for branch_id, name in enumerate(branch_result["names"]):
        result = pd.concat(
            [
                result,
                pd.DataFrame.from_dict(
                    {
                        "name": [name] * len(branch_result["time"]),
                        "time": np.array(branch_result["time"]),
                        "flow_in": np.array(
                            branch_result["flow_in"][branch_id]
                        ),
                        "flow_out": np.array(
                            branch_result["flow_out"][branch_id]
                        ),
                        "pressure_in": np.array(
                            branch_result["pressure_in"][branch_id]
                        ),
                        "pressure_out": np.array(
                            branch_result["pressure_out"][branch_id]
                        ),
                    }
                ),
            ]
        )
    return result


def run_from_file(input_file, output_file):
    """Run the svZeroDSolver from file.

    Args:
        input_file: Input file with configuration.
        output_file: Output file with configuration.
    """
    with open(input_file) as ff:
        config = json.load(ff)
    result = run_from_config(config)
    with open(output_file, "w") as ff:
        json.dump(result, ff)


@click.command()
@click.option(
    "--input_file",
    help="Path to the svZeroDSolver input file.",
    required=True,
    type=str,
)
@click.option(
    "--output_file",
    help="Path to the svZeroDSolver output file.",
    required=True,
    type=str,
)
def _run_from_from_sys_args(input_file, output_file):
    """Run the svZeroDSolver."""
    run_from_file(input_file, output_file)
