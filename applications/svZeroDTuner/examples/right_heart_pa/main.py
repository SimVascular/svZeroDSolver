"""
sv0D Tuning Framework - right_heart_pa example.

This script provides three modes:
1. BASELINE MODE: Run the initial model and save all results for inspection
2. SENSITIVITY MODE: Run correlation-based sensitivity screening
3. OPTIMIZE MODE: Run optimization using targets specified in tuning yaml

Recommended command-line workflow:
    Baseline: python -c 'from main import run_baseline; run_baseline("model.json")'
    Optimize: svzerodtuner optimize tuning_differential_evolution.yaml
"""

import os
import sys
import numpy as np
import pandas as pd
import pysvzerod

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

from svzerodtuner.sv0d_tuner import SV0DTuner
from svzerodtuner.visualization import plot_simulation_results
from svzerodtuner.sensitivity import SensitivityAnalyzer


def run_baseline(config_file):
    """
    Run the baseline simulation and save all results for user inspection.
    """
    print("="*70)
    print("BASELINE SIMULATION")
    print("="*70)
    print()

    if not os.path.exists(config_file):
        print(f"ERROR: Config file not found: {config_file}")
        return

    print(f"Running simulation with: {config_file}")
    print()

    try:
        solver = pysvzerod.Solver(config_file)
        solver.run()
        print("\u2713 Simulation completed successfully\n")
    except Exception as e:
        print(f"ERROR running simulation: {e}")
        return

    times = solver.get_times()
    full_results = solver.get_full_result()
    result_names = full_results['name'].unique()

    print(f"Found {len(result_names)} output variables")
    print()

    results_data = {'time': times}
    summary_stats = []

    for name in result_names:
        try:
            values = solver.get_single_result(name)
            results_data[name] = values
            stats = {
                'output_name': name,
                'min': np.min(values),
                'max': np.max(values),
                'mean': np.mean(values),
                'std': np.std(values)
            }
            summary_stats.append(stats)
        except Exception as e:
            print(f"Warning: Could not extract {name}: {e}")

    output_dir = 'baseline_results'
    os.makedirs(output_dir, exist_ok=True)

    results_df = pd.DataFrame(results_data)
    baseline_file = os.path.join(output_dir, 'baseline_results.csv')
    results_df.to_csv(baseline_file, index=False)
    print(f"\u2713 Saved full time series to: {baseline_file}")

    summary_df = pd.DataFrame(summary_stats)
    summary_file = os.path.join(output_dir, 'baseline_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f"\u2713 Saved summary statistics to: {summary_file}")
    print()

    print("="*70)
    print("BASELINE RESULTS SUMMARY")
    print("="*70)
    print()
    print(f"{'Output Variable':<40} {'Min':>12} {'Max':>12} {'Mean':>12}")
    print("-"*70)

    for stats in summary_stats:
        print(f"{stats['output_name']:<40} {stats['min']:>12.4e} "
              f"{stats['max']:>12.4e} {stats['mean']:>12.4e}")

    print()
    print("Generating plots...")
    plot_simulation_results(results_df, summary_df, output_dir, title_prefix="Baseline")

    print()
    print("="*70)
    print("NEXT STEPS:")
    print("="*70)
    print(f"1. Inspect {output_dir}/baseline_results.csv and baseline_summary.csv")
    print(f"2. View plots in {output_dir}/ to visualize the outputs")
    print("3. Choose which outputs you want to target")
    print("4. Update tuning yaml with your desired targets")
    print("5. Run svzerodtuner optimize <tuning_config.yaml>")
    print("="*70)
    print()


def run_optimization(config_file):
    """
    Run optimization using targets specified in config_file.
    """
    print("="*70)
    print("OPTIMIZATION")
    print("="*70)
    print()

    if not os.path.exists(config_file):
        print(f"ERROR: Config file not found: {config_file}")
        print(f"Please create {config_file} with your optimization settings.")
        return

    print(f"Using configuration: {config_file}")
    print()

    try:
        tuner = SV0DTuner(config_file)
    except Exception as e:
        print(f"ERROR loading configuration: {e}")
        return

    print("Starting optimization...")
    print("="*70)
    print()

    try:
        tuner.optimize()
    except Exception as e:
        print(f"ERROR during optimization: {e}")
        return


def run_sensitivity(config_file):
    """
    Run correlation-based sensitivity screening.
    """
    print("="*70)
    print("SENSITIVITY ANALYSIS")
    print("="*70)
    print()

    if not os.path.exists(config_file):
        print(f"ERROR: Config file not found: {config_file}")
        print(f"Please create {config_file} with your sensitivity analysis settings.")
        return

    print(f"Using configuration: {config_file}")
    print()

    try:
        analyzer = SensitivityAnalyzer(config_file)
    except Exception as e:
        print(f"ERROR loading configuration: {e}")
        return

    print("Starting sensitivity analysis...")
    print("="*70)
    print()

    try:
        results = analyzer.run()
    except Exception as e:
        print(f"ERROR during sensitivity analysis: {e}")
        return

    try:
        analyzer.save_results()
    except Exception as e:
        print(f"ERROR saving results: {e}")
        return


def main():
    # Choose one of the following:
    # baseline_config = "model.json"
    # optimize_config = "tuning_differential_evolution.yaml"
    # optimize_config = "tuning_nelder_mead.yaml"

    # Uncomment the mode you want to run:

    # run_baseline("model.json")
    run_optimization("tuning_differential_evolution.yaml")
    # run_optimization("tuning_nelder_mead.yaml")
    # run_sensitivity("sensitivity.yaml")
    pass


if __name__ == "__main__":
    main()
