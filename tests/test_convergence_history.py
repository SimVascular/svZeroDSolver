"""Tests for the Newton convergence history reported for every time step."""

import json
import os
import subprocess
import sys

import pandas as pd
import pytest

import svzerod

this_file_dir = os.path.abspath(os.path.dirname(__file__))

# Number of steps the steady initial condition takes (see Solver::setup_initial)
NUM_STEADY_STEPS = 31

# Test case with a non-trivial number of Newton iterations per time step
CASE = "valve_tanh"


def read_case(name=CASE, **parameters):
    """Read a test case and override its simulation parameters."""
    with open(os.path.join(this_file_dir, "cases", name + ".json")) as ff:
        config = json.load(ff)

    config["simulation_parameters"].update(parameters)

    return config


def run_case(name=CASE, **parameters):
    """Run a test case and return the solver and its configuration."""
    config = read_case(name, **parameters)
    solver = svzerod.Solver(config)
    solver.run()

    return solver, config


def test_history_covers_every_step():
    """All quantities are recorded for every step of the whole run."""
    solver, config = run_case()

    times = solver.get_step_times()
    iterations = solver.get_nonlinear_iterations()
    residuals = solver.get_step_residuals()
    equations = solver.get_step_residual_equations()
    phases = solver.get_step_phases()

    # steady initialization plus all steps of the transient simulation
    num_cycles = config["simulation_parameters"]["number_of_cardiac_cycles"]
    num_per_cycle = config["simulation_parameters"][
        "number_of_time_pts_per_cardiac_cycle"
    ]
    num_steps = NUM_STEADY_STEPS + (num_per_cycle - 1) * num_cycles

    assert len(times) == num_steps
    for quantity in [iterations, residuals, equations, phases]:
        assert len(quantity) == num_steps

    # the steady initialization comes first
    assert phases[:NUM_STEADY_STEPS] == ["steady_initial"] * NUM_STEADY_STEPS
    assert set(phases[NUM_STEADY_STEPS:]) == {"transient"}

    # every step needed at least one Newton update and stayed below the limit
    max_iter = config["simulation_parameters"].get(
        "maximum_nonlinear_iterations", 30
    )
    assert min(iterations) >= 1
    assert max(iterations) <= max_iter

    # this case is not solved in a single iteration everywhere
    assert max(iterations) > 1

    # every step converged and reports an equation of the model
    assert max(residuals) < config["simulation_parameters"]["absolute_tolerance"]
    assert "unknown" not in equations


def test_history_without_steady_initial():
    """Without a steady initial condition, only transient steps are recorded."""
    solver, config = run_case(steady_initial=False)

    phases = solver.get_step_phases()

    assert set(phases) == {"transient"}
    assert len(phases) == len(solver.get_step_times())

    num_cycles = config["simulation_parameters"]["number_of_cardiac_cycles"]
    num_per_cycle = config["simulation_parameters"][
        "number_of_time_pts_per_cardiac_cycle"
    ]
    assert len(phases) == (num_per_cycle - 1) * num_cycles


def test_history_as_dataframe():
    """The history is also available as a table."""
    solver, _ = run_case()

    history = solver.get_convergence_history()

    assert list(history.columns) == [
        "time",
        "phase",
        "iterations",
        "residual",
        "equation",
    ]
    assert len(history) == len(solver.get_step_times())
    assert history["iterations"].tolist() == solver.get_nonlinear_iterations()
    assert history["equation"].tolist() == solver.get_step_residual_equations()


def test_residual_belongs_to_the_returned_solution():
    """A step converging with its last update is not reported as a failure.

    The Newton loop checks the residual before each update, so the update of
    the last iteration is only checked after the loop. Running with exactly as
    many iterations as the hardest step of the reference run needs must still
    converge everywhere.
    """
    solver, config = run_case()
    atol = config["simulation_parameters"]["absolute_tolerance"]
    max_iter = max(solver.get_nonlinear_iterations())

    limited, _ = run_case(
        maximum_nonlinear_iterations=max_iter, max_iter_error_to_warning=True
    )

    assert max(limited.get_nonlinear_iterations()) == max_iter
    assert max(limited.get_step_residuals()) < atol
    assert limited.get_nonlinear_iterations() == solver.get_nonlinear_iterations()


def test_too_few_iterations_are_reported():
    """A step that really does not converge is still reported."""
    config = read_case(maximum_nonlinear_iterations=1)

    with pytest.raises(RuntimeError, match="non-linear iterations"):
        solver = svzerod.Solver(config)
        solver.run()


def run_cli(config, tmp_path):
    """Run a configuration the way it is run from the terminal."""
    input_file = os.path.join(tmp_path, "input.json")
    output_file = os.path.join(tmp_path, "output.csv")

    with open(input_file, "w") as ff:
        json.dump(config, ff)

    executable = os.path.join(this_file_dir, "..", "Release", "svzerodsolver")
    if os.name == "nt":
        executable += ".exe"

    if os.path.exists(executable):
        subprocess.run([executable, input_file, output_file], check=True)
    else:
        # the Python module exposes the same entry point
        argv = sys.argv
        sys.argv = ["svzerodsolver", input_file, output_file]
        try:
            svzerod.run_simulation_cli()
        finally:
            sys.argv = argv

    return output_file


def test_history_written_from_the_command_line(tmp_path):
    """The history of a run started from the terminal ends up in a file."""
    config = read_case(output_convergence_history=True)
    output_file = run_cli(config, tmp_path)

    history_file = os.path.join(tmp_path, "output_convergence.csv")
    assert os.path.exists(history_file)

    history = pd.read_csv(history_file)
    reference = svzerod.Solver(config)
    reference.run()

    assert list(history.columns) == [
        "time",
        "phase",
        "iterations",
        "residual",
        "equation",
    ]
    assert history["iterations"].tolist() == reference.get_nonlinear_iterations()
    assert history["phase"].tolist() == reference.get_step_phases()


def test_history_not_written_unless_requested(tmp_path):
    """The extra file is only written when the configuration asks for it."""
    run_cli(read_case(), tmp_path)

    assert not os.path.exists(os.path.join(tmp_path, "output_convergence.csv"))
