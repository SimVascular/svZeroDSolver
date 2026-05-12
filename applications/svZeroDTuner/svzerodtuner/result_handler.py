"""
Result Handler for sv0D Tuning Framework

Handles storing optimization history, saving results, and generating reports.
"""

import json
import os
import numpy as np
import pandas as pd
from typing import List, Dict, Optional
from .visualization import (
    plot_objective_history,
    plot_parameter_evolution,
    plot_target_comparison,
    plot_simulation_results,
    create_optimization_report
)


# Tolerance for "close to bound" warning: within this fraction of the range from a bound
_BOUND_TOLERANCE = 0.01


def _check_params_near_bounds(
    params: Dict[str, float], parameters: List[Dict]
) -> List[tuple]:
    """
    Check which parameters are close to their bounds.
    
    Returns:
        List of (param_name, value, bound_type, bound_value) for params near bounds
    """
    near_bounds = []
    for p in parameters:
        name = p['name']
        if name not in params:
            continue
        bounds = p.get('bounds')
        if not bounds or len(bounds) != 2:
            continue
        lower, upper = float(bounds[0]), float(bounds[1])
        value = float(params[name])
        range_val = upper - lower
        if range_val <= 0:
            continue
        tol = _BOUND_TOLERANCE * range_val
        if value <= lower + tol:
            near_bounds.append((name, value, 'min', lower))
        elif value >= upper - tol:
            near_bounds.append((name, value, 'max', upper))
    return near_bounds


class ResultHandler:
    """
    Handles storage and reporting of optimization results.
    """
    
    def __init__(self, output_dir: str, save_history: bool = True, save_plots: bool = True):
        """
        Initialize result handler.
        
        Args:
            output_dir: Output directory for results
            save_history: Whether to save optimization history
            save_plots: Whether to generate plots
        """
        self.output_dir = output_dir
        self.save_history = save_history
        self.save_plots = save_plots
        
        os.makedirs(output_dir, exist_ok=True)
    
    def save_history_json(self, history: List[Dict], filename: str = "history.json"):
        """
        Save optimization history to JSON file.
        
        Args:
            history: Optimization history list
            filename: Output filename
        """
        if not self.save_history or not history:
            if not history:
                print("Note: No optimization history to save")
            return
        
        # Save history in optimization_history subfolder
        history_dir = os.path.join(self.output_dir, 'optimization_history')
        os.makedirs(history_dir, exist_ok=True)
        filepath = os.path.join(history_dir, filename)
        
        # Convert numpy arrays to lists for JSON serialization
        history_serializable = []
        for entry in history:
            entry_copy = entry.copy()
            if 'parameters' in entry_copy:
                entry_copy['parameters'] = {
                    k: float(v) if isinstance(v, (np.number, np.ndarray)) else v
                    for k, v in entry_copy['parameters'].items()
                }
            history_serializable.append(entry_copy)
        
        with open(filepath, 'w') as f:
            json.dump(history_serializable, f, indent=2)
        
        print(f"Saved optimization history to {filepath}")
    
    def save_history_csv(self, history: List[Dict], filename: str = "history.csv"):
        """
        Save optimization history to CSV file.
        
        Args:
            history: Optimization history list
            filename: Output filename
        """
        if not self.save_history or not history:
            return
        
        # Save history in optimization_history subfolder
        history_dir = os.path.join(self.output_dir, 'optimization_history')
        os.makedirs(history_dir, exist_ok=True)
        filepath = os.path.join(history_dir, filename)
        
        rows = []
        for entry in history:
            step = entry.get('iteration', len(rows))
            row = {
                'iteration': step, 
                'objective': entry['objective']
            }
            row.update(entry['parameters'])
            rows.append(row)
        
        df = pd.DataFrame(rows)
        df.to_csv(filepath, index=False)
        
        print(f"Saved optimization history CSV to {filepath}")
    
    def save_final_parameters(
        self,
        param_names: List[str],
        param_values: np.ndarray,
        filename: str = "final_parameters.json"
    ):
        """
        Save final optimized parameters.
        
        Args:
            param_names: Parameter names
            param_values: Parameter values
            filename: Output filename
        """
        filepath = os.path.join(self.output_dir, filename)
        
        params_dict = {
            name: float(value) for name, value in zip(param_names, param_values)
        }
        
        with open(filepath, 'w') as f:
            json.dump(params_dict, f, indent=2)
        
        print(f"Saved final parameters to {filepath}")
    
    def save_final_config(
        self,
        config_dict: Dict,
        filename: str = "optimized_config.json"
    ):
        """
        Save final optimized sv0D configuration.
        
        Args:
            config_dict: sv0D configuration dictionary
            filename: Output filename
        """
        filepath = os.path.join(self.output_dir, filename)
        
        with open(filepath, 'w') as f:
            json.dump(config_dict, f, indent=2)
        
        print(f"Saved optimized configuration to {filepath}")
    
    def generate_plots(
        self,
        history: List[Dict],
        param_names: List[str],
        targets: Optional[List[Dict]] = None,
        simulated_values: Optional[Dict] = None,
        optimized_results_df: Optional[pd.DataFrame] = None,
        optimized_summary_df: Optional[pd.DataFrame] = None
    ):
        """
        Generate optimization visualization plots.
        
        Args:
            history: Optimization history
            param_names: Parameter names
            targets: Target specifications (optional)
            simulated_values: Simulated values dictionary (optional)
            optimized_results_df: DataFrame with optimized simulation time series (optional)
            optimized_summary_df: DataFrame with optimized simulation statistics (optional)
        """
        if not self.save_plots:
            return
        
        # Create subfolder for optimization history/evolution plots
        history_dir = os.path.join(self.output_dir, 'optimization_history')
        os.makedirs(history_dir, exist_ok=True)
        
        if history:
            # Plot objective history
            plot_objective_history(
                history,
                output_file=os.path.join(history_dir, 'objective_history.png')
            )
            
            # Plot parameter evolution
            if param_names:
                plot_parameter_evolution(
                    history,
                    param_names,
                    output_file=os.path.join(history_dir, 'parameter_evolution.png')
                )
        else:
            print("Note: Skipping history plots (no history available)")
        
        # Plot target comparison if available
        if targets and simulated_values:
            plot_target_comparison(
                targets,
                simulated_values,
                output_file=os.path.join(self.output_dir, 'target_comparison.png')
            )
        
        # Plot all simulation results if DataFrames are provided
        if optimized_results_df is not None and optimized_summary_df is not None:
            # Create subfolder for simulation results
            sim_output_dir = os.path.join(self.output_dir, 'optimized_simulation')
            os.makedirs(sim_output_dir, exist_ok=True)
            
            # Save optimized results CSVs
            optimized_results_df.to_csv(os.path.join(sim_output_dir, 'optimized_results.csv'), index=False)
            optimized_summary_df.to_csv(os.path.join(sim_output_dir, 'optimized_summary.csv'), index=False)
            print(f"Saved optimized results to: {sim_output_dir}/")
            
            # Plot optimized results
            plot_simulation_results(
                optimized_results_df, 
                optimized_summary_df, 
                sim_output_dir,
                title_prefix="Optimized"
            )
    
    def create_summary_report(
        self,
        history: List[Dict],
        param_names: List[str],
        best_value: float,
        best_params: np.ndarray,
        parameters: Optional[List[Dict]] = None,
        targets: Optional[List[Dict]] = None,
        simulated_values: Optional[Dict] = None,
        optimized_results_df: Optional[pd.DataFrame] = None,
        optimized_summary_df: Optional[pd.DataFrame] = None,
        timing_info: Optional[Dict] = None
    ):
        """
        Create comprehensive summary report.
        
        Args:
            history: Optimization history
            param_names: Parameter names
            best_value: Best objective value
            best_params: Best parameter values
            parameters: Parameter definitions with bounds (optional, for bound warnings)
            targets: Target specifications (optional)
            simulated_values: Simulated values dictionary (optional)
            optimized_results_df: DataFrame with optimized simulation time series (optional)
            optimized_summary_df: DataFrame with optimized simulation statistics (optional)
            timing_info: Dictionary with timing statistics (optional)
        """
        # Save history
        self.save_history_json(history)
        self.save_history_csv(history)
        
        # Save final parameters
        self.save_final_parameters(param_names, best_params)
        
        # Save timing information if provided
        if timing_info:
            history_dir = os.path.join(self.output_dir, 'optimization_history')
            os.makedirs(history_dir, exist_ok=True)
            timing_file = os.path.join(history_dir, 'timing_info.json')
            with open(timing_file, 'w') as f:
                json.dump(timing_info, f, indent=2)
            print(f"Saved timing information to {timing_file}")
        
        # Generate plots
        self.generate_plots(
            history,
            param_names,
            targets=targets,
            simulated_values=simulated_values,
            optimized_results_df=optimized_results_df,
            optimized_summary_df=optimized_summary_df
        )
        
        # Print summary
        print("\n" + "="*60)
        print("OPTIMIZATION SUMMARY")
        print("="*60)
        if timing_info and 'n_iterations' in timing_info:
            print(f"Total iterations: {timing_info['n_iterations']}")
            n_ev = timing_info.get('n_evaluations')
            print(f"Total function evaluations: {n_ev if n_ev is not None else 'not available (interrupted)'}")
        else:
            print(f"Total function evaluations: {len(history)}")
        print(f"Best objective value: {best_value:.6e}")
        print(f"\nBest parameters:")
        for name, value in zip(param_names, best_params):
            print(f"  {name:<50} {value:.6e}")
        # Warn if any best parameter is close to its bounds
        if parameters:
            near_bounds = _check_params_near_bounds(
                dict(zip(param_names, best_params)), parameters
            )
            if near_bounds:
                print("\nWARNING: The following best parameters are near their bounds:")
                for name, value, bound_type, bound_value in near_bounds:
                    print(f"  {name}={value:.6e} is near {bound_type} bound ({bound_value:.6e})")
        print("="*60)
