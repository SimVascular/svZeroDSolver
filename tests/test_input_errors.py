"""Check that invalid input configurations fail with messages that name the
offending block, config key, and value, and list the valid alternatives.

Each test loads a working case from tests/cases/, introduces one mistake, and
checks the RuntimeError raised while the model is loaded.
"""

import copy
import json
import os
import re

import pytest

import svzerod

this_file_dir = os.path.abspath(os.path.dirname(__file__))

# Coupled case with an inlet coupling block connected to the ClosedLoopRCR
# boundary condition "RCR_aorta" and an outlet coupling block connected to
# "ClosedLoopHeartAndPulmonary".
COUPLED_CASE = "coupledBlock_closedLoopHeart_singleVessel.json"


def load_case(name):
    with open(os.path.join(this_file_dir, "cases", name)) as ff:
        return json.load(ff)


@pytest.fixture(scope="module")
def coupled_config():
    return load_case(COUPLED_CASE)


def get_coupling_block(config, name):
    return next(
        block
        for block in config["external_solver_coupling_blocks"]
        if block["name"] == name
    )


def assert_load_error(config, message):
    with pytest.raises(RuntimeError, match=re.escape(message)):
        svzerod.simulate(config)


def test_misspelled_bc_type_connected_to_coupling_block(coupled_config):
    # The scenario from GitHub issue #209: previously this failed with a
    # generic "connection type for inlet external_coupling_block is invalid".
    config = copy.deepcopy(coupled_config)
    config["boundary_conditions"][0]["bc_type"] = "ClosedLoopRCRR"
    assert_load_error(
        config,
        "Boundary condition 'RCR_aorta' has invalid bc_type 'ClosedLoopRCRR'. "
        "Valid bc_type values are: ",
    )
    # The list of valid types comes from the block factory; spot-check it.
    with pytest.raises(RuntimeError) as err:
        svzerod.simulate(config)
    valid = str(err.value).split("Valid bc_type values are: ")[1]
    for bc_type in ["ClosedLoopRCR", "FLOW", "IMPEDANCE", "RCR", "RESISTANCE"]:
        assert bc_type in valid
    for not_a_bc_type in ["BloodVessel", "ClosedLoopHeartAndPulmonary"]:
        assert not_a_bc_type not in valid


def test_misspelled_bc_type_without_coupling():
    config = load_case("pulsatileFlow_R_RCR.json")
    bc = next(b for b in config["boundary_conditions"] if b["bc_type"] == "RCR")
    bc["bc_type"] = "RCRR"
    assert_load_error(
        config,
        f"Boundary condition '{bc['bc_name']}' has invalid bc_type 'RCRR'.",
    )


def test_coupling_block_unknown_connected_block(coupled_config):
    config = copy.deepcopy(coupled_config)
    get_coupling_block(config, "external_inlet")["connected_block"] = "RCR_aorto"
    assert_load_error(
        config,
        "External coupling block 'external_inlet' has connected_block "
        "'RCR_aorto', which is not the bc_name of a boundary condition",
    )


def test_coupling_block_unsupported_connected_type(coupled_config):
    config = copy.deepcopy(coupled_config)
    get_coupling_block(config, "external_inlet")[
        "connected_block"
    ] = "ClosedLoopHeartAndPulmonary"
    assert_load_error(
        config,
        "External coupling block 'external_inlet' (location 'inlet') has "
        "connected_block 'ClosedLoopHeartAndPulmonary' of type "
        "'ClosedLoopHeartAndPulmonary', which is not supported at this "
        "location. Supported types for location 'inlet' are: RESISTANCE, RCR, "
        "IMPEDANCE, ClosedLoopRCR, CORONARY, ClosedLoopCoronaryLeft, "
        "ClosedLoopCoronaryRight, BloodVessel.",
    )


def test_coupling_block_invalid_location(coupled_config):
    config = copy.deepcopy(coupled_config)
    get_coupling_block(config, "external_inlet")["location"] = "inelt"
    assert_load_error(
        config,
        "External coupling block 'external_inlet' has invalid location "
        "'inelt'. Valid location values are: inlet, outlet.",
    )
