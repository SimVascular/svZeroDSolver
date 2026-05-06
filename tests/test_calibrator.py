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


def test_pulsatile_R_only_recovers_from_y_and_t():
    """Calibrate only ``R_poiseuille`` from a pulsatile forward simulation
    using only ``y`` and ``t``; ``dy`` is reconstructed internally from ``y``
    via the generalized-alpha relation. The fixture was produced by
    ``pysvzerod.simulate`` with R=100, C=1e-4, L=1; the calibrator should
    recover R close to 100 (the reconstruction is first-order accurate so the
    error is small but non-zero)."""
    testfile = os.path.join(
        this_file_dir,
        "cases",
        "vmr",
        "input",
        "pulsatileFlow_R_calibrate_R_only.json",
    )
    result, _ = execute_pysvzerod(testfile, "calibrator")
    R = result["vessels"][0]["zero_d_element_values"]["R_poiseuille"]
    assert np.isclose(R, 100.0, rtol=1e-2), f"R = {R}, expected ~100"


def test_calibration_rejects_time_dependent_block_without_t(tmp_path):
    """The calibrator must refuse a model containing a time-dependent block
    when the input lacks an observation time vector ``t``."""
    with open(
        os.path.join(this_file_dir, "cases", "steadyFlow_calibration.json")
    ) as ff:
        cfg = json.load(ff)
    cfg["vessels"][0]["zero_d_element_type"] = "ChamberSphere"
    cfg["vessels"][0]["calibrate"] = ["W1"]

    testfile = tmp_path / "rejects_time_dependent.json"
    with open(testfile, "w") as ff:
        json.dump(cfg, ff)

    with pytest.raises(RuntimeError, match="time-dependent"):
        execute_pysvzerod(str(testfile), "calibrator")


@pytest.mark.parametrize(
    "test_case",
    [
        "0080_0001_calibrate_from_0d",
        "0104_0001_calibrate_from_0d",
        "0140_2001_calibrate_from_0d",
        # Calibrates only R_poiseuille via per-block ``calibrate`` fields,
        # starting from the reference values for C, L and
        # stenosis_coefficient and zeroed R_poiseuille. The calibrator should
        # recover R_poiseuille while leaving the other parameters untouched.
        "0104_0001_calibrate_R_only",
    ],
)
def test_calibration_vmr(test_case):
    """Calibrate a model from the vascular model repository and check that
    every parameter matches the corresponding reference."""
    test = os.path.join(
        this_file_dir, "cases", "vmr", "input", f"{test_case}.json"
    )
    model_id = test_case[:9]
    reference_file = os.path.join(
        this_file_dir, "cases", "vmr", "reference", f"{model_id}_optimal_from_0d.json"
    )
    with open(reference_file) as ff:
        reference = json.load(ff)

    result, _ = execute_pysvzerod(test, "calibrator")

    for i, vessel in enumerate(reference["vessels"]):
        for key, value in vessel["zero_d_element_values"].items():
            assert np.isclose(
                result["vessels"][i]["zero_d_element_values"][key],
                value,
                rtol=RTOL_PRES,
            )

    for i, junction in enumerate(reference["junctions"]):
        if "junction_values" in junction:
            for key, value in junction["junction_values"].items():
                assert np.allclose(
                    result["junctions"][i]["junction_values"][key],
                    value,
                    rtol=RTOL_PRES,
                )
