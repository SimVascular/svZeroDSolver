import copy
import json
import os
import pytest

import numpy as np

from .utils import execute_pysvzerod, RTOL_PRES

this_file_dir = os.path.abspath(os.path.dirname(__file__))


def test_steady_flow_calibration():
    testfile = os.path.join(this_file_dir, "cases", "steadyFlow_calibration.json")

    result, _ = execute_pysvzerod(testfile, "calibrator")

    calibrated_parameters = result["vessels"][0]["zero_d_element_values"]

    assert np.isclose(
        np.mean(calibrated_parameters["R_poiseuille"]), 100, rtol=RTOL_PRES
    )
    assert np.isclose(np.mean(calibrated_parameters["C"]), 0.0001, rtol=RTOL_PRES)
    assert np.isclose(np.mean(calibrated_parameters["L"]), 1.0, rtol=RTOL_PRES)
    assert np.isclose(
        np.mean(calibrated_parameters["stenosis_coefficient"]), 0.0, rtol=RTOL_PRES
    )


def _load_vmr_case(model_id):
    """Load the VMR input (with y/dy) and reference (with optimal params)."""
    with open(
        os.path.join(
            this_file_dir, "cases", "vmr", "input", f"{model_id}_calibrate_from_0d.json"
        )
    ) as ff:
        input_cfg = json.load(ff)
    with open(
        os.path.join(
            this_file_dir,
            "cases",
            "vmr",
            "reference",
            f"{model_id}_optimal_from_0d.json",
        )
    ) as ff:
        reference = json.load(ff)
    return input_cfg, reference


def _build_R_only_config(input_cfg, reference):
    """Compose a calibrator config that starts from the reference values, with
    every R_poiseuille zeroed; the y/dy/calibration_parameters come from the
    forward-simulation input file."""
    cfg = copy.deepcopy(reference)
    cfg["y"] = input_cfg["y"]
    cfg["dy"] = input_cfg["dy"]
    cfg["calibration_parameters"] = copy.deepcopy(input_cfg["calibration_parameters"])

    for vessel in cfg["vessels"]:
        vessel["zero_d_element_values"]["R_poiseuille"] = 0.0
    for junction in cfg["junctions"]:
        if "junction_values" in junction and "R_poiseuille" in junction["junction_values"]:
            junction["junction_values"]["R_poiseuille"] = [
                0.0 for _ in junction["junction_values"]["R_poiseuille"]
            ]
    return cfg


def _assert_recovers_reference(result, reference):
    for ref_vessel, res_vessel in zip(reference["vessels"], result["vessels"]):
        for key, ref_value in ref_vessel["zero_d_element_values"].items():
            res_value = res_vessel["zero_d_element_values"][key]
            assert np.isclose(res_value, ref_value, atol=1e-9, rtol=RTOL_PRES), (
                f"Vessel {ref_vessel['vessel_name']} parameter {key}: "
                f"expected {ref_value}, got {res_value}"
            )

    for ref_junction, res_junction in zip(reference["junctions"], result["junctions"]):
        if "junction_values" not in ref_junction:
            continue
        for key, ref_values in ref_junction["junction_values"].items():
            res_values = res_junction["junction_values"][key]
            assert np.allclose(res_values, ref_values, atol=1e-9, rtol=RTOL_PRES), (
                f"Junction {ref_junction['junction_name']} parameter {key}: "
                f"expected {ref_values}, got {res_values}"
            )


def test_calibration_vmr_R_only_global(tmp_path):
    """Calibrate only R_poiseuille via the global ``calibrate`` field.

    Starts from the calibrated reference (so C, L, stenosis_coefficient are at
    the ground truth), zeros every R_poiseuille and asks the calibrator to
    recover R only via ``calibration_parameters.calibrate=["R_poiseuille"]``.
    """
    model_id = "0104_0001"
    input_cfg, reference = _load_vmr_case(model_id)
    cfg = _build_R_only_config(input_cfg, reference)
    cfg["calibration_parameters"]["calibrate"] = ["R_poiseuille"]

    testfile = tmp_path / f"{model_id}_calibrate_R_only_global.json"
    with open(testfile, "w") as ff:
        json.dump(cfg, ff)

    result, _ = execute_pysvzerod(str(testfile), "calibrator")
    _assert_recovers_reference(result, reference)


def test_calibration_vmr_R_only_per_block(tmp_path):
    """Calibrate only R_poiseuille via per-block ``calibrate`` fields, with no
    global default. Every vessel and junction independently lists the
    parameters it wants calibrated."""
    model_id = "0104_0001"
    input_cfg, reference = _load_vmr_case(model_id)
    cfg = _build_R_only_config(input_cfg, reference)

    for vessel in cfg["vessels"]:
        vessel["calibrate"] = ["R_poiseuille"]
    for junction in cfg["junctions"]:
        if "junction_values" in junction:
            junction["calibrate"] = ["R_poiseuille"]

    testfile = tmp_path / f"{model_id}_calibrate_R_only_per_block.json"
    with open(testfile, "w") as ff:
        json.dump(cfg, ff)

    result, _ = execute_pysvzerod(str(testfile), "calibrator")
    _assert_recovers_reference(result, reference)


def test_calibration_vmr_block_overrides_global(tmp_path):
    """Per-block ``calibrate`` must override the global default. Sets a
    nonsense global default of ``["C"]`` (which would calibrate the wrong
    parameter and fail the recovery test) and overrides every block to
    ``["R_poiseuille"]``. The global must be ignored where overridden."""
    model_id = "0104_0001"
    input_cfg, reference = _load_vmr_case(model_id)
    cfg = _build_R_only_config(input_cfg, reference)
    cfg["calibration_parameters"]["calibrate"] = ["C"]

    for vessel in cfg["vessels"]:
        vessel["calibrate"] = ["R_poiseuille"]
    for junction in cfg["junctions"]:
        if "junction_values" in junction:
            junction["calibrate"] = ["R_poiseuille"]

    testfile = tmp_path / f"{model_id}_calibrate_block_overrides.json"
    with open(testfile, "w") as ff:
        json.dump(cfg, ff)

    result, _ = execute_pysvzerod(str(testfile), "calibrator")
    _assert_recovers_reference(result, reference)


@pytest.mark.parametrize("model_id", ["0080_0001", "0104_0001", "0140_2001"])
def test_calibration_vmr(model_id):
    """Test actual models from the vascular model repository."""
    with open(
        os.path.join(
            this_file_dir, "cases", "vmr", "input", f"{model_id}_calibrate_from_0d.json"
        )
    ) as ff:
        reference = json.load(ff)

    test = os.path.join(
        this_file_dir, "cases", "vmr", "input", f"{model_id}_calibrate_from_0d.json"
    )

    result, _ = execute_pysvzerod(test, "calibrator")

    for i, vessel in enumerate(reference["vessels"]):
        for key, value in vessel["zero_d_element_values"].items():
            np.isclose(
                result["vessels"][i]["zero_d_element_values"][key],
                value,
                rtol=RTOL_PRES,
            )

    for i, junction in enumerate(reference["junctions"]):
        if "junction_values" in junction:
            for key, value in junction["junction_values"].items():
                np.allclose(
                    result["junctions"][i]["junction_values"][key],
                    value,
                    rtol=RTOL_PRES,
                )
