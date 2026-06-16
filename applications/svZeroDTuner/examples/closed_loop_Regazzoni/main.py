"""
sv0D Tuning Framework - Regazzoni closed-loop model example.

This script provides three modes:
1. BASELINE MODE: Run the initial model and save all results for inspection
2. SENSITIVITY MODE: Run correlation-based sensitivity screening
3. OPTIMIZE MODE: Run optimization using targets specified in tuning.yaml

Recommended command-line workflow:
    Baseline: python -c 'from main import run_baseline; run_baseline("model.json")'
    Optimize: svzerodtuner optimize tuning_nelder_mead.yaml
    Sensitivity: svzerodtuner sensitivity sensitivity.yaml
"""

import os
import sys
import numpy as np
import pandas as pd
import pysvzerod

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import json
from svzerodtuner.sv0d_tuner import SV0DTuner
from svzerodtuner.visualization import plot_simulation_results
from svzerodtuner.sensitivity import SensitivityAnalyzer


def run_baseline(config_file):
    """
    Run the baseline simulation and save all results for user inspection.
    
    This function:
    - Runs the initial model.json simulation
    - Saves all available outputs to baseline_results.csv
    - Displays summary statistics (min, max, mean) for each output
    - User can then inspect these results and specify targets in tuning_config.yaml
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
    
    # Run baseline simulation
    try:
        solver = pysvzerod.Solver(config_file)
        solver.run()
        print("✓ Simulation completed successfully\n")
    except Exception as e:
        print(f"ERROR running simulation: {e}")
        return
    
    # Extract all results
    times = solver.get_times()
    full_results = solver.get_full_result()
    result_names = full_results['name'].unique()
    
    print(f"Found {len(result_names)} output variables")
    print()
    
    # Create results DataFrame with all outputs
    results_data = {'time': times}
    summary_stats = []
    
    for name in result_names:
        try:
            values = solver.get_single_result(name)
            results_data[name] = values
            
            # Calculate statistics
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
    
    # Create baseline_results directory
    output_dir = 'baseline_results'
    os.makedirs(output_dir, exist_ok=True)
    
    # Save full results to CSV
    results_df = pd.DataFrame(results_data)
    baseline_file = os.path.join(output_dir, 'baseline_results.csv')
    results_df.to_csv(baseline_file, index=False)
    print(f"✓ Saved full time series to: {baseline_file}")
    
    # Save summary statistics
    summary_df = pd.DataFrame(summary_stats)
    summary_file = os.path.join(output_dir, 'baseline_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f"✓ Saved summary statistics to: {summary_file}")
    print()
    
    # Display summary statistics
    print("="*70)
    print("BASELINE RESULTS SUMMARY")
    print("="*70)
    print()
    print(f"{'Output Variable':<40} {'Min':>12} {'Max':>12} {'Mean':>12}")
    print("-"*70)
    
    for stats in summary_stats:
        print(f"{stats['output_name']:<40} {stats['min']:>12.4e} "
              f"{stats['max']:>12.4e} {stats['mean']:>12.4e}")
    
    # Generate plots
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
    print("4. Update tuning.yaml with your desired targets")
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
    
    # Initialize tuner
    try:
        tuner = SV0DTuner(config_file)
    except Exception as e:
        print(f"ERROR loading configuration: {e}")
        return
    
    # Run optimization
    print("Starting optimization...")
    print("="*70)
    print()
    
    try:
        results = tuner.optimize()
    except Exception as e:
        print(f"ERROR during optimization: {e}")
        return

def run_sensitivity(config_file):
    """
    Run correlation-based sensitivity screening.
    
    This function:
    - Performs global sensitivity analysis on specified parameters
    - Computes first-order and total-order screening scores
    - Identifies which parameters most influence each quantity of interest
    - Saves results, plots, and summary statistics
    """
    print("="*70)
    print("SENSITIVITY ANALYSIS")
    print("="*70)
    print()
    
    if not os.path.exists(config_file):
        print(f"ERROR: Config file not found: {config_file}")
        print(f"Please create {config_file} with your sensitivity analysis settings.")
        print(f"See sensitivity.yaml.example for template.")
        return
    
    print(f"Using configuration: {config_file}")
    print()
    
    # Initialize sensitivity analyzer
    try:
        analyzer = SensitivityAnalyzer(config_file)
    except Exception as e:
        print(f"ERROR loading configuration: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Run sensitivity analysis
    try:
        results = analyzer.run()
    except Exception as e:
        print(f"ERROR during sensitivity analysis: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Save results
    try:
        analyzer.save_results()
    except Exception as e:
        print(f"ERROR saving results: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Print summary
    print()
    print("="*70)
    print("SENSITIVITY ANALYSIS COMPLETE")
    print("="*70)
    print()
    print("Summary of Results:")
    print("-"*70)
    
    for qoi_key, qoi_results in results.items():
        print(f"\n{qoi_key}:")
        print(f"  Range: [{qoi_results['min']:.4e}, {qoi_results['max']:.4e}]")
        print(f"  Mean ± Std: {qoi_results['mean']:.4e} ± {qoi_results['std']:.4e}")
        print(f"\n  Most influential parameters (first-order indices):")
        
        # Sort by influence
        first_order = qoi_results['first_order']
        sorted_params = sorted(first_order.items(), key=lambda x: abs(x[1]), reverse=True)
        
        for param, value in sorted_params:
            print(f"    {param:<30} {value:>8.4f}")
    
    print()
    print(f"Results saved to: {analyzer.output_config.get('directory', 'sensitivity_results')}")
    print("="*70)
    print()


def main():
    """
    Main function.
    
    INSTRUCTIONS:
    ============
    Uncomment ONE of the following modes to run:
    
    MODE 1: BASELINE    - Run initial simulation and save results for inspection
    MODE 2: OPTIMIZE    - Run optimization with targets from tuning.yaml
    MODE 3: SENSITIVITY - Run global sensitivity analysis with sensitivity.yaml
    """
    # Change to script directory
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # ============================================================================
    # SELECT MODE: Uncomment ONE of the following
    # ============================================================================
    
    #run_baseline("model.json")                     # MODE 1: Run baseline and save results
    #run_sensitivity("sensitivity.yaml")            # MODE 2: Run sensitivity analysis
    run_optimization("tuning_nelder_mead.yaml")     # MODE 3: Run optimization with tuning.yaml
    
    
    # ============================================================================


if __name__ == "__main__":
    main()
