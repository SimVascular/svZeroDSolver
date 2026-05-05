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


def test_calibration_vmr_R_only(tmp_path):
    """Calibrate only R_poiseuille on a VMR model while holding C, L and
    stenosis_coefficient constant at their ground-truth values.

    Starts from the calibrated reference (so C, L, stenosis_coefficient are at
    the ground truth), zeros out every R_poiseuille in vessels and junctions,
    and runs the calibrator with ``calibrate=["R_poiseuille"]``. The y/dy
    observations come from the corresponding input file. The calibrator should
    recover the reference R_poiseuille values and leave every other parameter
    untouched.
    """
    model_id = "0104_0001"

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

    cfg = copy.deepcopy(reference)
    cfg["y"] = input_cfg["y"]
    cfg["dy"] = input_cfg["dy"]
    cfg["calibration_parameters"] = copy.deepcopy(input_cfg["calibration_parameters"])
    cfg["calibration_parameters"]["calibrate"] = ["R_poiseuille"]

    # Zero every R_poiseuille; leave C, L, stenosis_coefficient at ground truth
    for vessel in cfg["vessels"]:
        vessel["zero_d_element_values"]["R_poiseuille"] = 0.0
    for junction in cfg["junctions"]:
        if "junction_values" in junction and "R_poiseuille" in junction["junction_values"]:
            junction["junction_values"]["R_poiseuille"] = [
                0.0 for _ in junction["junction_values"]["R_poiseuille"]
            ]

    testfile = tmp_path / f"{model_id}_calibrate_R_only.json"
    with open(testfile, "w") as ff:
        json.dump(cfg, ff)

    result, _ = execute_pysvzerod(str(testfile), "calibrator")

    # R_poiseuille is recovered; other parameters are unchanged
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
