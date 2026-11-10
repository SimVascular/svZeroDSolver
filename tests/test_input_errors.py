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
