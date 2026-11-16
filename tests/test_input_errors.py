"""Check that invalid input configurations fail with messages that name the
offending block, config key, and value, and list the valid alternatives.

Each test loads a working case from tests/cases/, introduces one mistake, and
checks the error raised while the model is loaded.
"""

import copy
import json
import os
import subprocess
from tempfile import TemporaryDirectory

import pytest

import svzerod

this_file_dir = os.path.abspath(os.path.dirname(__file__))

# Coupled case with an inlet coupling block "external_inlet" connected to the
# ClosedLoopRCR boundary condition "RCR_aorta", and an outlet coupling block
# "external_outlet" connected to "ClosedLoopHeartAndPulmonary".
COUPLED_CASE = "coupledBlock_closedLoopHeart_singleVessel"


def load_case(name):
    with open(os.path.join(this_file_dir, "cases", name + ".json")) as ff:
        return json.load(ff)


def get_coupling_block(config, name):
    return next(
        block
        for block in config["external_solver_coupling_blocks"]
        if block["name"] == name
    )


def load_error_message(config):
    """Run `config`, which must fail to load, and return the error message.

    With --coverage, run the coverage-instrumented Release executable (like
    utils.execute_svzerod) so these error paths count toward the coverage
    report; the solver then reports the error on stderr and exits with 1.
    """
    if pytest.coverage:
        with TemporaryDirectory() as tempdir:
            config_file = os.path.join(tempdir, "config.json")
            with open(config_file, "w") as ff:
                json.dump(config, ff)
            exe = os.path.join(this_file_dir, "..", "Release", "svzerodsolver")
            run = subprocess.run(
                [exe, config_file, os.path.join(tempdir, "out.csv")],
                capture_output=True,
                text=True,
            )
        assert run.returncode == 1, run.stdout + run.stderr
        return run.stderr
    with pytest.raises(RuntimeError) as err:
        svzerod.simulate(config)
    return str(err.value)


def test_misspelled_bc_type_connected_to_coupling_block():
    # The scenario from GitHub issue #209: previously this failed with a
    # generic "connection type for inlet external_coupling_block is invalid".
    config = load_case(COUPLED_CASE)
    config["boundary_conditions"][0]["bc_type"] = "ClosedLoopRCRR"
    message = load_error_message(config)
    assert (
        "Block 'RCR_aorta' in boundary_conditions has invalid bc_type "
        "'ClosedLoopRCRR'. Valid bc_type values are: " in message
    )
    # The list of valid types comes from the block factory; spot-check it.
    # Compare whole names: "RCR" is a substring of "ClosedLoopRCR".
    listed = message.split("Valid bc_type values are: ")[1]
    valid = set(listed.split(".")[0].split(", "))
    for bc_type in ["ClosedLoopRCR", "FLOW", "IMPEDANCE", "RCR", "RESISTANCE"]:
        assert bc_type in valid
    for not_a_bc_type in ["BloodVessel", "ClosedLoopHeartAndPulmonary"]:
        assert not_a_bc_type not in valid


@pytest.mark.parametrize(
    "case, component, index, key, valid",
    [
        ("pulsatileFlow_R_RCR", "boundary_conditions", 1, "bc_type", None),
        (
            "steadyFlow_R_steadyPressure",
            "vessels",
            0,
            "zero_d_element_type",
            "BloodVessel, BloodVesselCRL, BloodVesselRC",
        ),
        (
            "steadyFlow_confluenceR_R",
            "junctions",
            0,
            "junction_type",
            "BloodVesselJunction, NORMAL_JUNCTION, resistive_junction",
        ),
        (
            "closedLoopHeart_singleVessel",
            "closed_loop_blocks",
            0,
            "closed_loop_type",
            "ClosedLoopHeartAndPulmonary",
        ),
        ("chamber_sphere", "valves", 0, "type", "PiecewiseValve, ValveTanh"),
        ("chamber_sphere", "chambers", 0, "type", None),
        (COUPLED_CASE, "external_solver_coupling_blocks", 0, "type", None),
    ],
)
def test_misspelled_block_type(case, component, index, key, valid):
    """Each config section names the block and key of a misspelled type."""
    config = load_case(case)
    block = config[component][index]
    block[key] += "_typo"
    name = next(
        (block[k] for k in ["bc_name", "vessel_name", "junction_name", "name"]
         if k in block),
        f"entry {index}",
    )
    message = load_error_message(config)
    assert (
        f"Block '{name}' in {component} has invalid {key} '{block[key]}'. "
        f"Valid {key} values are: " in message
    )
    if valid is not None:
        assert f"Valid {key} values are: {valid}." in message
    if component == "external_solver_coupling_blocks":
        # Every bc_type loads as a coupling type, but only FLOW and PRESSURE
        # exchange data with the external solver, so the error says so.
        assert message.rstrip().endswith(
            ". External solvers normally couple through FLOW or PRESSURE."
        )


def test_coupling_block_unknown_connected_block():
    config = load_case(COUPLED_CASE)
    get_coupling_block(config, "external_inlet")["connected_block"] = "RCR_aorto"
    assert (
        "External coupling block 'external_inlet' has connected_block "
        "'RCR_aorto', which is not the bc_name of a boundary condition"
        in load_error_message(config)
    )


def test_coupling_block_unsupported_connected_type_inlet():
    config = load_case(COUPLED_CASE)
    get_coupling_block(config, "external_inlet")[
        "connected_block"
    ] = "ClosedLoopHeartAndPulmonary"
    assert (
        "External coupling block 'external_inlet' (location 'inlet') has "
        "connected_block 'ClosedLoopHeartAndPulmonary' of type "
        "'ClosedLoopHeartAndPulmonary', which is not supported at this "
        "location. Supported types for location 'inlet' are: RESISTANCE, RCR, "
        "IMPEDANCE, ClosedLoopRCR, CORONARY, ClosedLoopCoronaryLeft, "
        "ClosedLoopCoronaryRight, BloodVessel." in load_error_message(config)
    )


def test_coupling_block_unsupported_connected_type_outlet():
    config = load_case(COUPLED_CASE)
    config["boundary_conditions"][0]["bc_type"] = "RESISTANCE"
    get_coupling_block(config, "external_outlet")["connected_block"] = "RCR_aorta"
    assert (
        "External coupling block 'external_outlet' (location 'outlet') has "
        "connected_block 'RCR_aorta' of type 'RESISTANCE', which is not "
        "supported at this location. Supported types for location 'outlet' "
        "are: ClosedLoopRCR, ClosedLoopHeartAndPulmonary, BloodVessel."
        in load_error_message(config)
    )


def test_coupling_block_invalid_location():
    config = load_case(COUPLED_CASE)
    get_coupling_block(config, "external_inlet")["location"] = "inelt"
    assert (
        "External coupling block 'external_inlet' has invalid location "
        "'inelt'. Valid location values are: inlet, outlet."
        in load_error_message(config)
    )


# Bifurcation with vessels 0, 1, 2, junction "J0" from vessel 0 to vessels 1
# and 2, and the boundary conditions INFLOW, OUT1 and OUT2.
BIFURCATION_CASE = "steadyFlow_bifurcationR_R1"


def add_unconnected_vessel(config, vessel_id, vessel_name):
    """Append a copy of vessel "branch1_seg0" with no boundary conditions."""
    vessel = copy.deepcopy(config["vessels"][1])
    vessel["vessel_id"] = vessel_id
    vessel["vessel_name"] = vessel_name
    vessel.pop("boundary_conditions", None)
    config["vessels"].append(vessel)
    return vessel


def test_vessel_without_any_connection():
    """A vessel left out of every junction used to segfault the solver.

    Blocks index their nodes' degrees of freedom by position, so a vessel with
    no connections read past the end of its own variable list.
    """
    config = load_case(BIFURCATION_CASE)
    add_unconnected_vessel(config, 3, "orphan")
    message = load_error_message(config)
    assert "Invalid connectivity in the model:" in message
    assert (
        "Block 'orphan' (BloodVessel) has 0 inlet connections and 0 outlet "
        "connections, but needs exactly 1 inlet and exactly 1 outlet."
        in message
    )


def test_vessel_with_only_one_connection():
    """A vessel connected on one side names the side that is missing."""
    config = load_case(BIFURCATION_CASE)
    config["junctions"][0]["outlet_vessels"] = [1]
    assert (
        "Block 'branch2_seg0' (BloodVessel) has 0 inlet connections and 1 "
        "outlet connection, but needs exactly 1 inlet and exactly 1 outlet."
        in load_error_message(config)
    )


def test_junction_without_inlet():
    """A junction keeps its own requirement of at least one inlet."""
    config = load_case(BIFURCATION_CASE)
    config["junctions"][0]["inlet_vessels"] = []
    config["vessels"][0].pop("boundary_conditions")
    config["boundary_conditions"] = config["boundary_conditions"][1:]
    assert (
        "Block 'J0' (NORMAL_JUNCTION) has 0 inlet connections and 2 outlet "
        "connections, but needs at least 1 inlet and at least 1 outlet."
        in load_error_message(config)
    )


def test_boundary_condition_without_connection():
    """A boundary condition no vessel refers to is reported, not ignored."""
    config = load_case(BIFURCATION_CASE)
    config["boundary_conditions"].append(
        {"bc_name": "UNUSED", "bc_type": "RESISTANCE",
         "bc_values": {"Pd": 100.0, "R": 100.0}}
    )
    assert (
        "Block 'UNUSED' (RESISTANCE) has 0 connections, but a boundary "
        "condition needs exactly 1, as the inlet or the outlet of one other "
        "block." in load_error_message(config)
    )


def test_junction_with_unknown_vessel_id():
    """A vessel_id no vessel defines used to leave the junction dangling."""
    config = load_case(BIFURCATION_CASE)
    config["junctions"][0]["outlet_vessels"] = [1, 5]
    assert (
        "Junction 'J0' lists vessel_id 5 in outlet_vessels, but no vessel has "
        "that vessel_id. Defined vessel_id values are: 0, 1, 2."
        in load_error_message(config)
    )


def test_junction_without_connection_keys():
    """Half a connection pair leaves a junction with no connections at all."""
    config = load_case(BIFURCATION_CASE)
    config["junctions"][0].pop("outlet_vessels")
    assert (
        "Junction 'J0' does not define its connections. Provide either both "
        "inlet_vessels and outlet_vessels, or both inlet_blocks and "
        "outlet_blocks. Keys found: inlet_vessels." in load_error_message(config)
    )


def test_junction_lists_only_the_first_vessel_ids():
    """Models have hundreds of vessels, so the list of ids is cut short."""
    config = load_case(BIFURCATION_CASE)
    for vessel_id in range(3, 15):
        add_unconnected_vessel(config, vessel_id, f"extra{vessel_id}")
    config["junctions"][0]["inlet_vessels"] = [99]
    assert (
        "Junction 'J0' lists vessel_id 99 in inlet_vessels, but no vessel has "
        "that vessel_id. Defined vessel_id values are: 0, 1, 2, 3, 4, 5, 6, 7, "
        "8, 9, and 5 more." in load_error_message(config)
    )


def load_calibrator_error_message(config):
    """Run `config` through the calibrator, which must fail, and return the
    error message.

    The calibrator reads the `vessels` and `junctions` sections itself, so it
    needs its own check that a junction's vessel_id error reaches the user.
    """
    if pytest.coverage:
        with TemporaryDirectory() as tempdir:
            config_file = os.path.join(tempdir, "config.json")
            with open(config_file, "w") as ff:
                json.dump(config, ff)
            exe = os.path.join(this_file_dir, "..", "Release",
                               "svzerodcalibrator")
            run = subprocess.run(
                [exe, config_file, os.path.join(tempdir, "out.json")],
                capture_output=True,
                text=True,
            )
        assert run.returncode == 1, run.stdout + run.stderr
        return run.stderr
    with pytest.raises(RuntimeError) as err:
        svzerod.calibrate(config)
    return str(err.value)


def test_calibrator_junction_with_unknown_vessel_id():
    """A vessel_id no vessel defines used to reach the calibrator as a block
    lookup for the empty name."""
    config = load_case("vmr/input/0080_0001_calibrate_from_0d")
    config["junctions"][0]["outlet_vessels"] = [1, 9999]
    assert (
        "Junction 'J0' lists vessel_id 9999 in outlet_vessels, but no vessel "
        "has that vessel_id. Defined vessel_id values are: "
        in load_calibrator_error_message(config)
    )


# The non-coverage CI job runs pytest before building Release/, so tests that
# need an executable only run where one exists.
CALIBRATOR_EXE = os.path.join(
    this_file_dir, "..", "Release", "svzerodcalibrator"
)


@pytest.mark.skipif(
    not os.path.exists(CALIBRATOR_EXE),
    reason="needs the svzerodcalibrator executable",
)
def test_calibrator_malformed_json_file():
    """A file that isn't valid JSON used to terminate the calibrator.

    Only the executable reads a file; the Python binding takes a dict.
    """
    with TemporaryDirectory() as tempdir:
        config_file = os.path.join(tempdir, "config.json")
        with open(config_file, "w") as ff:
            ff.write('{"vessels": [,]}')
        run = subprocess.run(
            [CALIBRATOR_EXE, config_file, os.path.join(tempdir, "out.json")],
            capture_output=True,
            text=True,
        )
    assert run.returncode == 1, run.stdout + run.stderr
    assert (
        f"[svzerodcalibrator] Error: Parsing the input file '{config_file}' "
        "has failed." in run.stderr
    )
