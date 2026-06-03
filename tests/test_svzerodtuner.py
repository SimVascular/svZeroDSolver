import importlib
import json
import sys
import types
import tomllib
from pathlib import Path

import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_pyproject_declares_svzerodtuner_test_dependencies():
    with open(REPO_ROOT / "pyproject.toml", "rb") as f:
        pyproject = tomllib.load(f)

    installed_deps = {
        dep.split(";")[0].split("[")[0].strip().lower()
        for dep in pyproject["project"]["dependencies"]
    }
    installed_deps.update(
        dep.split(";")[0].split("[")[0].strip().lower()
        for dep in pyproject["dependency-groups"]["dev"]
    )

    assert {"numpy", "pandas", "scipy", "pyyaml", "matplotlib"} <= installed_deps


def test_parameter_handler_rejects_unknown_parameter_names(tmp_path):
    model_path = tmp_path / "model.json"
    model_path.write_text(
        json.dumps(
            {
                "chambers": [{"name": "LV", "values": {"Emax": 2.0}}],
                "vessels": [
                    {
                        "vessel_name": "AORTA",
                        "zero_d_element_values": {"R_poiseuille": 1.0},
                        "vessel_length": 10.0,
                    }
                ],
                "valves": [{"name": "MV", "params": {"Rmin": 0.1}}],
                "boundary_conditions": [{"bc_name": "OUT", "bc_values": {"R": 5.0}}],
                "simulation_parameters": {"density": 1.06},
                "initial_condition": {"pressure": 80.0},
            }
        )
    )

    from applications.svZeroDTuner.svzerodtuner.parameter_handler import ParameterHandler

    handler = ParameterHandler(str(model_path))
    handler.set_parameter("LV.Emax", 3.0)
    assert handler.get_parameter("LV.Emax") == 3.0

    with pytest.raises(ValueError, match="Parameter 'LV.Unknown' not found"):
        handler.set_parameter("LV.Unknown", 1.0)

    with pytest.raises(ValueError, match="Parameter 'AORTA.fake_field' not found"):
        handler.set_parameter("AORTA.fake_field", 1.0)


def test_sensitivity_analyzer_validates_parameter_names_up_front(tmp_path):
    sys.modules.setdefault("pysvzerod", types.SimpleNamespace(Solver=object))
    sensitivity_module = importlib.import_module(
        "applications.svZeroDTuner.svzerodtuner.sensitivity"
    )

    model_path = tmp_path / "model.json"
    model_path.write_text(json.dumps({"chambers": [{"name": "LV", "values": {"Emax": 2.0}}]}))

    config_path = tmp_path / "sensitivity.yaml"
    config_path.write_text(
        "\n".join(
            [
                "model:",
                f"  config_file: {model_path}",
                "parameters:",
                "  - name: LV.Unknown",
                "    bounds: [1.0, 3.0]",
                "quantities_of_interest:",
                "  - name: lv_pressure_max",
                "    expression: np.max(pressure:LV)",
            ]
        )
    )

    with pytest.raises(
        ValueError, match="Sensitivity parameter 'LV.Unknown' not found"
    ):
        sensitivity_module.SensitivityAnalyzer(str(config_path))


def test_negative_relative_bounds_are_ordered_for_scalar_targets():
    from applications.svZeroDTuner.svzerodtuner.objective import ObjectiveFunction

    objective = ObjectiveFunction(
        targets=[
            {
                "name": "venous_pressure",
                "type": "scalar",
                "target_value": -100.0,
                "relative_bounds": "10%",
            }
        ],
        norm="L1",
    )

    target = objective.targets[0]
    assert float(target["range_lo"][0]) == pytest.approx(-110.0)
    assert float(target["range_hi"][0]) == pytest.approx(-90.0)
    assert objective.compute({"venous_pressure": -100.0}) == pytest.approx(0.0)


def test_config_handler_rejects_reversed_time_series_target_range(tmp_path):
    model_path = tmp_path / "model.json"
    model_path.write_text(json.dumps({"chambers": [{"name": "LV", "values": {"Emax": 2.0}}]}))

    target_path = tmp_path / "target.csv"
    target_path.write_text("time,value\n0.0,1.0\n1.0,2.0\n")

    config_path = tmp_path / "tuning.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "model": {"config_file": str(model_path)},
                "parameters": [{"name": "LV.Emax", "bounds": [1.0, 3.0]}],
                "targets": [
                    {
                        "name": "lv_volume",
                        "type": "time_series",
                        "expression": "Vc:LV",
                        "target_file": str(target_path),
                        "target_range": [5.0, -5.0],
                    }
                ],
                "objective": {"norm": "L1"},
                "optimization": {"algorithm": "Nelder-Mead"},
            }
        )
    )

    from applications.svZeroDTuner.svzerodtuner.config_handler import ConfigHandler

    with pytest.raises(ValueError, match="target_range must have min < max"):
        ConfigHandler(str(config_path))


def test_sensitivity_results_use_screening_labels_and_filenames(tmp_path, monkeypatch):
    sys.modules.setdefault("pysvzerod", types.SimpleNamespace(Solver=object))
    sensitivity_module = importlib.import_module(
        "applications.svZeroDTuner.svzerodtuner.sensitivity"
    )

    analyzer = object.__new__(sensitivity_module.SensitivityAnalyzer)
    analyzer.output_config = {}
    analyzer.sample_data = [{"sample_id": 0, "LV.Emax": 2.0, "qoi": 1.5}]
    analyzer.parameters = [{"name": "LV.Emax"}]
    analyzer.n_samples = 1
    analyzer.results = {
        "qoi": {
            "analysis_type": "correlation_screening",
            "sampler": "sobol_sequence",
            "first_order_metric": "squared_pearson_correlation",
            "total_order_metric": "binned_conditional_variance_screening",
            "first_order": {"LV.Emax": 0.25},
            "total_order": {"LV.Emax": 0.5},
            "mean": 1.5,
            "std": 0.0,
            "min": 1.5,
            "max": 1.5,
        }
    }
    monkeypatch.setattr(analyzer, "_create_visualizations", lambda output_path: None)

    analyzer.save_results(str(tmp_path))

    assert (tmp_path / "screening_indices.json").exists()
    assert not (tmp_path / "sobol_indices.json").exists()
    summary = (tmp_path / "summary.txt").read_text()
    assert "correlation-based sensitivity screening" in summary
    assert "Sobol indices" not in summary
    assert "First-order screening score" in summary


def test_sensitivity_run_resets_sample_data_each_time(monkeypatch):
    sys.modules.setdefault("pysvzerod", types.SimpleNamespace(Solver=object))
    sensitivity_module = importlib.import_module(
        "applications.svZeroDTuner.svzerodtuner.sensitivity"
    )

    analyzer = object.__new__(sensitivity_module.SensitivityAnalyzer)
    analyzer.parameters = [{"name": "LV.Emax", "bounds": [1.0, 2.0]}]
    analyzer.quantities_of_interest = [{"name": "qoi"}]
    analyzer.sensitivity_config = {"n_samples": 1}
    analyzer.output_config = {}
    analyzer.n_samples = 1
    analyzer.results = {"stale": {}}
    analyzer.sample_data = [{"sample_id": 999, "LV.Emax": -1.0, "qoi": -1.0}]

    samples = [
        {"qoi": 1.0},
        {"qoi": 2.0},
    ]

    def fake_get_qoi(_param_values):
        return samples.pop(0)

    monkeypatch.setattr(analyzer, "_get_simulated_quantities_of_interest", fake_get_qoi)

    first_run = analyzer.run()
    assert len(analyzer.sample_data) == 1
    assert analyzer.sample_data[0]["sample_id"] == 0
    assert analyzer.sample_data[0]["qoi"] == pytest.approx(1.0)
    assert "stale" not in first_run

    samples.append({"qoi": 2.0})
    second_run = analyzer.run()
    assert len(analyzer.sample_data) == 1
    assert analyzer.sample_data[0]["sample_id"] == 0
    assert analyzer.sample_data[0]["qoi"] == pytest.approx(2.0)
    assert "stale" not in second_run
