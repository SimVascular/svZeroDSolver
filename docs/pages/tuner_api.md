@page tuner_api svZeroDTuner API Reference

[TOC]

# CLI Commands

Command-line entrypoint:

```bash
svzerodtuner <command> <config.yaml>
```

Supported commands:

- `optimize <config>`: run optimization from tuning YAML
- `run <config>`: alias for `optimize`
- `sensitivity-analysis <config>`: run sensitivity analysis from YAML
- `sensitivity <config>`: alias for `sensitivity-analysis`

# Python Entry Points

Primary optimization API:

- `SV0DTuner(config_file: str)`
- `SV0DTuner.optimize() -> Dict`
- `SV0DTuner.evaluate(param_values: Optional[Dict[str, float]] = None) -> Dict`
- `run_optimization(config_file: str) -> Dict`

Primary sensitivity API:

- `SensitivityAnalyzer(config_file: str)`
- `SensitivityAnalyzer.run() -> Dict`
- `SensitivityAnalyzer.save_results(output_dir: Optional[str] = None) -> None`
- `run_sensitivity_analysis(config_file: str) -> Dict`

Supporting components:

- `ConfigHandler` (YAML loading/validation)
- `ParameterHandler` (model parameter get/set by name)
- `OutputExtractor` (output extraction by exact solver output name)
- `OptimizerWrapper` (SciPy optimizer orchestration)
- `ObjectiveFunction` (range-based target penalties)
- `Expression` (expression parsing/evaluation for targets and QoIs)

# Data Contracts

Optimization return dictionary includes:

- `success` (`bool`)
- `message` (`str`)
- `best_value` (`float`)
- `best_params` (`Dict[str, float]`)
- `history` (`List[Dict]`)
- `result` (SciPy `OptimizeResult`)
- `interrupted` (`bool`)

Evaluation return dictionary includes:

- `objective_value` (`float`)
- `simulated_values` (`Dict`)
- `parameters` (`Dict[str, float]`)

Sensitivity result dictionary is keyed by QoI name and contains:

- `first_order` (`Dict[str, float]`)
- `total_order` (`Dict[str, float]`)
- `mean`, `std`, `min`, `max` (`float`)

# Source Pointers

- CLI: `applications/svZeroDTuner/svzerodtuner/cli.py`
- Optimization orchestrator: `applications/svZeroDTuner/svzerodtuner/sv0d_tuner.py`
- Sensitivity: `applications/svZeroDTuner/svzerodtuner/sensitivity.py`
- Configuration: `applications/svZeroDTuner/svzerodtuner/config_handler.py`
- Objective and optimizer internals:
  - `applications/svZeroDTuner/svzerodtuner/objective.py`
  - `applications/svZeroDTuner/svzerodtuner/optimizer.py`
