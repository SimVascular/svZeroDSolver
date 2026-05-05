import json
import os

import numpy as np
import pytest

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


@pytest.mark.parametrize(
    "case",
    [
        # Selects parameters via the global ``calibration_parameters.calibrate``.
        "0104_0001_calibrate_R_only_global",
        # Selects parameters via per-block ``calibrate`` fields, no global default.
        "0104_0001_calibrate_R_only_per_block",
        # Per-block ``calibrate`` overrides a misleading global default.
        "0104_0001_calibrate_R_only_block_overrides",
    ],
)
def test_calibration_R_only(case):
    """Calibrate only ``R_poiseuille`` on a VMR model while holding C, L and
    stenosis_coefficient at their ground-truth values. The test fixtures start
    from the calibrated reference (so non-R parameters are at the optimum) with
    every R_poiseuille zeroed; the calibrator should recover R_poiseuille and
    leave the rest untouched. The reference values to compare against live in
    ``tests/cases/vmr/reference/0104_0001_optimal_from_0d.json``.
    """
    testfile = os.path.join(this_file_dir, "cases", "vmr", "input", f"{case}.json")
    reference_file = os.path.join(
        this_file_dir, "cases", "vmr", "reference", "0104_0001_optimal_from_0d.json"
    )
    with open(reference_file) as ff:
        reference = json.load(ff)

    result, _ = execute_pysvzerod(testfile, "calibrator")

    for ref_vessel, res_vessel in zip(reference["vessels"], result["vessels"]):
        for key, ref_value in ref_vessel["zero_d_element_values"].items():
            res_value = res_vessel["zero_d_element_values"][key]
            assert np.isclose(res_value, ref_value, atol=1e-9, rtol=RTOL_PRES), (
                f"Vessel {ref_vessel['vessel_name']} parameter {key}: "
                f"expected {ref_value}, got {res_value}"
            )

    for ref_junction, res_junction in zip(
        reference["junctions"], result["junctions"]
    ):
        if "junction_values" not in ref_junction:
            continue
        for key, ref_values in ref_junction["junction_values"].items():
            res_values = res_junction["junction_values"][key]
            assert np.allclose(res_values, ref_values, atol=1e-9, rtol=RTOL_PRES), (
                f"Junction {ref_junction['junction_name']} parameter {key}: "
                f"expected {ref_values}, got {res_values}"
            )
