"""This module holds the main execution routines of svZeroDSolver."""
import json

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

    zero_d_results_branch = format_results_to_dict(
        time_steps, y_list, block_list
    )
    return zero_d_results_branch


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
