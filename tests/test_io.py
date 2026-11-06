import numpy as np
import pytest

import svzerod

from .utils import run_test_case_by_name, RTOL_FLOW, RTOL_PRES, get_test_case_by_name

# use coarse absolute tolerances for gradient calculation
# we're comparing gradients from gen-alpha in svZeroDSolver with central differences in np.gradient
ATOL = {"f": 0.6, "p": 800.0}
ATOL_MEAN = {"f": 0.01, "p": 20.0}


def test_pulsatile_flow_r_rcr_mean():
    # time-dependent results
    res_time = run_test_case_by_name("pulsatileFlow_R_RCR")

    # time-averaged results
    res_mean = run_test_case_by_name("pulsatileFlow_R_RCR_mean")

    # compare all time-dependent results to their average
    tol = {"f": RTOL_FLOW, "p": RTOL_PRES}
    for k, res in res_time.items():
        assert np.isclose(np.mean(res[0]), res_mean[k][0][0], rtol=tol[k[0]])


def test_pulsatile_flow_r_rcr_mean_variable():
    # time-dependent results
    res_time = run_test_case_by_name(
        "pulsatileFlow_R_RCR_variable", output_variable_based=True
    )

    # time-averaged results
    res_mean = run_test_case_by_name(
        "pulsatileFlow_R_RCR_mean_variable", output_variable_based=True
    )

    # fields and segments to test
    fields = ["flow", "pressure"]
    segments = {"in": "INFLOW:branch0_seg0", "out": "branch0_seg0:OUT"}

    tol = {"f": RTOL_FLOW, "p": RTOL_PRES}
    for f in fields:
        for k, v in segments.items():
            rt = np.mean(res_time[f + "_" + k])
            rm = res_mean["y"][res_mean["name"] == f + ":" + v]
            assert np.isclose(rt, rm, rtol=tol[f[0]])


def test_pulsatile_flow_r_rcr_derivative():
    results = run_test_case_by_name(
        "pulsatileFlow_R_RCR_derivative", output_variable_based=True
    )

    # time step
    dt = results["time"][1] - results["time"][0]

    # fields to test
    fields = ["pressure_in", "pressure_out", "flow_in", "flow_out"]
    for f in fields:
        ref = np.gradient(results[f], dt)
        res = results["d_" + f]
        assert np.all(np.isclose(ref, res, atol=ATOL[f[0]]))


def test_pulsatile_flow_r_rcr_derivative_variable():
    results = run_test_case_by_name(
        "pulsatileFlow_R_RCR_derivative_variable", output_variable_based=True
    )

    # time step
    dt = results["time"][1] - results["time"][0]

    # fields and segments to test
    fields = ["flow", "pressure"]
    segments = ["INFLOW:branch0_seg0", "branch0_seg0:OUT"]
    for f in fields:
        for v in segments:
            ids = results["name"] == f + ":" + v
            ref = np.gradient(results.y[ids], dt)
            res = results.ydot[ids]
            assert np.all(np.isclose(ref, res, atol=ATOL[f[0]]))


def test_pulsatile_flow_r_rcr_mean_derivative():
    # time-dependent results
    res_time = run_test_case_by_name("pulsatileFlow_R_RCR", output_variable_based=True)

    # time-averaged results
    res_mean = run_test_case_by_name(
        "pulsatileFlow_R_RCR_mean_derivative", output_variable_based=True
    )

    # time step
    dt = res_time["time"][1] - res_time["time"][0]

    # fields to test
    fields = ["pressure_in", "pressure_out", "flow_in", "flow_out"]
    for f in fields:
        ref = np.mean(np.gradient(res_time[f], dt))
        res = res_mean["d_" + f]
        assert np.isclose(ref, res, atol=ATOL_MEAN[f[0]])


def test_pulsatile_flow_r_rcr_mean_derivative_variable():
    # time-dependent results
    res_time = run_test_case_by_name(
        "pulsatileFlow_R_RCR_variable", output_variable_based=True
    )

    # time-averaged results
    res_mean = run_test_case_by_name(
        "pulsatileFlow_R_RCR_mean_derivative_variable", output_variable_based=True
    )

    # time step
    dt = res_time["time"][1] - res_time["time"][0]

    # fields and segments to test
    fields = ["flow", "pressure"]
    segments = {"in": "INFLOW:branch0_seg0", "out": "branch0_seg0:OUT"}
    for f in fields:
        for k, v in segments.items():
            rt = np.mean(np.gradient(res_time[f + "_" + k], dt))
            rm = res_mean["ydot"][res_mean["name"] == f + ":" + v]
            assert np.isclose(rt, rm, atol=ATOL_MEAN[f[0]])


def test_time_dependent_block():
    """FLOW and PRESSURE BCs given as an expression string follow it exactly."""
    res_flow = run_test_case_by_name("timeDep_Flow", output_variable_based=True)
    res_press = run_test_case_by_name("timeDep_Pressure", output_variable_based=True)

    name_flow = "flow:INLET:branch0_seg0"
    name_pres = "pressure:INLET:branch0_seg0"
    time = res_flow["time"][res_flow["name"] == name_flow].to_numpy()
    flow = res_flow["y"][res_flow["name"] == name_flow].to_numpy()
    pressure = res_press["y"][res_press["name"] == name_pres].to_numpy()

    # skip the initial transient
    expected = 2.0 * np.pi * np.cos(2.0 * np.pi * time)
    assert np.allclose(flow[15:], expected[15:], rtol=0.003)
    assert np.allclose(pressure[15:], expected[15:], rtol=0.003)


def run_expression_case(bc_values, case="timeDep_Flow", **sim_params):
    """Run an expression test case with the given inlet bc_values."""
    config = get_test_case_by_name(case)
    config["boundary_conditions"][0]["bc_values"] = bc_values
    config["simulation_parameters"].update(sim_params)
    return svzerod.simulate(config)


def test_expression_steady_initial_uses_cycle_mean():
    """The steady initial condition uses the expression's cycle average."""
    res = run_expression_case({"Q": "5.0 + sin(2 * pi * t)"}, steady_initial=True)
    pressure = res["y"][res["name"] == "pressure:INLET:branch0_seg0"].to_numpy()

    # steady state for the mean flow Q = 5 through R_poiseuille + R = 1100
    assert np.isclose(pressure[0], 5.0 * 1100.0, rtol=1e-6)


def test_expression_with_cardiac_period_from_simulation_parameters():
    """A model whose only time dependence is an expression can set cardiac_period."""
    res = run_expression_case(
        {"Q": "sin(4 * pi * t)"},
        cardiac_period=0.5,
        number_of_cardiac_cycles=2,
        number_of_time_pts_per_cardiac_cycle=51,
    )
    name = "flow:INLET:branch0_seg0"
    time = res["time"][res["name"] == name].to_numpy()
    flow = res["y"][res["name"] == name].to_numpy()
    assert np.isclose(time[-1], 1.0)

    # skip the initial transient
    assert np.allclose(flow[15:], np.sin(4 * np.pi * time[15:]), atol=5e-3)


@pytest.mark.parametrize(
    "case, bc_type, key",
    [("timeDep_Flow", "FLOW", "Q"), ("timeDep_Pressure", "PRESSURE", "P")],
)
@pytest.mark.parametrize(
    "bc_values, message",
    [
        ({"{key}": "sin(t"}, "Invalid expression in parameter {key} of {bc_type}"),
        ({"{key}": "x * t"}, "Invalid expression in parameter {key} of {bc_type}"),
        (
            {"{key}": "var x := 0; for (var i := 0; i < 10; i += 1) {{ x += t; }}; x"},
            "Invalid expression in parameter {key}",
        ),
        ({"{key}": "t := t * 2; t"}, "Invalid expression in parameter {key}"),
        ({"{key}": ""}, "Parameter {key} in {bc_type} block INLET is an empty string"),
        ({"{key}": "t", "t": [0.0, 1.0]}, "Parameter t in {bc_type} block INLET is unused"),
        ({}, "Array parameter {key} is mandatory in {bc_type} block INLET"),
        ({"{key}": [1.0, 1.0]}, "Array parameter t is mandatory in {bc_type} block INLET"),
    ],
)
def test_invalid_reference_bc_values(case, bc_type, key, bc_values, message):
    """Invalid FLOW and PRESSURE bc_values are rejected at load time."""
    bc_values = {
        k.format(key=key): v.format(key=key) if isinstance(v, str) else v
        for k, v in bc_values.items()
    }
    with pytest.raises(RuntimeError, match=message.format(key=key, bc_type=bc_type)):
        run_expression_case(bc_values, case=case)


def test_missing_mandatory_array_param():
    """Missing mandatory array parameter raises RuntimeError."""
    config = get_test_case_by_name("pulsatileFlow_R_coronary")
    del config["boundary_conditions"][1]["bc_values"]["Pim"]

    with pytest.raises(RuntimeError, match="Array parameter Pim is mandatory"):
        svzerod.simulate(config)


def test_missing_t_vector_for_array_param():
    """Mandatory array parameter present but missing t vector raises RuntimeError."""
    config = get_test_case_by_name("pulsatileFlow_R_coronary")
    del config["boundary_conditions"][1]["bc_values"]["t"]

    with pytest.raises(RuntimeError, match="Array parameter t is mandatory"):
        svzerod.simulate(config)
