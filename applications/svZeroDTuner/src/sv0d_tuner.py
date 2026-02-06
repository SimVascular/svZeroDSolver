"""
Main sv0D Tuning Framework

Orchestrates parameter optimization to match target output values.
"""

import os
import json
import time
import numpy as np
import pandas as pd
import pysvzerod
from typing import Dict, List, Optional, Callable

from .parameter_handler import ParameterHandler
from .output_extractor import OutputExtractor
from .objective import create_objective, ObjectiveFunction
from .optimizer import OptimizerWrapper
from .config_handler import ConfigHandler
from .result_handler import ResultHandler


class SV0DTuner:
    """
    Main framework class for tuning sv0D models.
    """
    
    def __init__(self, config_file: str):
        """
        Initialize sv0D tuner with configuration file.
        
        Args:
            config_file: Path to YAML configuration file
        """
        # Load configuration
        self.config_handler = ConfigHandler(config_file)
        self.config = self.config_handler.get_config()
        
        # Initialize components
        self.param_handler = ParameterHandler(
            self.config_handler.get_model_config_file()
        )
        
        # Get configuration sections
        self.parameters = self.config_handler.get_parameters()
        self.targets = self.config_handler.get_targets()
        self.objective_config = self.config_handler.get_objective_config()
        self.optimization_config = self.config_handler.get_optimization_config()
        self.output_config = self.config_handler.get_output_config()
        
        # Initialize optimizer
        self.optimizer = OptimizerWrapper(**self.optimization_config)
        
        # Initialize result handler
        self.result_handler = ResultHandler(
            output_dir=self.output_config['directory'],
            save_history=self.output_config.get('save_history', True),
            save_plots=self.output_config.get('save_plots', True)
        )
        
        # Create objective function
        self.objective_func = create_objective(
            objective_type=self.objective_config.get('type', 'weighted_l2'),
            targets=self.targets,
            normalize=self.objective_config.get('normalize', False),
            custom_function=self.objective_config.get('custom_function')
        )
        
        # State
        self.solver = None
        self.extractor = None
        self.history = []
        self.best_value = None
        self.best_params = None
    
    def _format_time(self, seconds: float) -> str:
        """
        Format time duration in human-readable format.
        
        Args:
            seconds: Time in seconds
            
        Returns:
            Formatted time string
        """
        if seconds < 1:
            return f"{seconds*1000:.1f} ms"
        elif seconds < 60:
            return f"{seconds:.2f} sec"
        elif seconds < 3600:
            minutes = int(seconds // 60)
            secs = seconds % 60
            return f"{minutes} min {secs:.1f} sec"
        else:
            hours = int(seconds // 3600)
            minutes = int((seconds % 3600) // 60)
            return f"{hours} hr {minutes} min"
    
    def _run_simulation(self, param_values: np.ndarray, return_full_results: bool = False):
        """
        Run sv0D simulation with given parameter values.
        
        Args:
            param_values: Array of parameter values
            return_full_results: If True, return tuple of (simulated_values, results_df, summary_df)
            
        Returns:
            If return_full_results=False: Dictionary of simulated output values
            If return_full_results=True: Tuple of (simulated_values dict, results_df, summary_df)
        """
        # Update parameters - ensure values are Python scalars
        param_names = [p['name'] for p in self.parameters]
        for name, value in zip(param_names, param_values):
            # Explicitly convert to Python float
            if isinstance(value, np.ndarray):
                value = float(value.item())
            elif not isinstance(value, (int, float)):
                value = float(value)
            self.param_handler.set_parameter(name, value)
        
        # Get updated config
        config_dict = self.param_handler.get_config()
        
        # Create and run solver directly with config dict
        # All numpy types have been converted by param_handler.get_config()
        self.solver = pysvzerod.Solver(config_dict)
        self.solver.run()
        
        # Create extractor
        self.extractor = OutputExtractor(self.solver)
        
        # Extract all target outputs
        simulated_values = {}
        times = self.extractor.get_times()
        
        # First, extract all unique outputs as time_series (needed for scalar extractions)
        unique_outputs = {}
        for target in self.targets:
            name = target['name']
            if name not in unique_outputs:
                # Always extract as time_series first (we can compute scalars from it)
                try:
                    time_series = self.extractor.extract(name, 'time_series')
                    unique_outputs[name] = {
                        'time_series': time_series,
                        'times': times
                    }
                except:
                    # If time_series extraction fails, try the requested type
                    unique_outputs[name] = {
                        'value': self.extractor.extract(name, target.get('type', 'time_series'))
                    }
        
        # Now process each target and compute the requested extraction type
        for target in self.targets:
            name = target['name']
            extraction_type = target.get('type', 'time_series')
            
            # Create unique key for targets with same name but different types
            target_key = f"{name}_{extraction_type}" if extraction_type != 'time_series' else name
            
            if name in unique_outputs and 'time_series' in unique_outputs[name]:
                # We have time series data, compute requested type
                ts_data = unique_outputs[name]['time_series']
                if extraction_type == 'time_series':
                    simulated_values[target_key] = ts_data
                    simulated_values[f'{target_key}_times'] = unique_outputs[name]['times']
                    # Also store under base name for backward compatibility
                    simulated_values[name] = ts_data
                    simulated_values[f'{name}_times'] = unique_outputs[name]['times']
                elif extraction_type == 'min':
                    simulated_values[target_key] = float(np.min(ts_data))
                    simulated_values[name] = simulated_values[target_key]  # Also store under base name
                elif extraction_type == 'max':
                    simulated_values[target_key] = float(np.max(ts_data))
                    simulated_values[name] = simulated_values[target_key]  # Also store under base name
                elif extraction_type == 'mean':
                    simulated_values[target_key] = float(np.mean(ts_data))
                    simulated_values[name] = simulated_values[target_key]  # Also store under base name
            else:
                # Use pre-extracted value
                simulated_values[target_key] = unique_outputs[name]['value']
        
        # If full results requested, extract all outputs
        if return_full_results:
            try:
                full_results = self.solver.get_full_result()
                result_names = full_results['name'].unique()
                
                # Create results DataFrame with all outputs
                results_data = {'time': times}
                summary_stats = []
                
                for name in result_names:
                    try:
                        values = self.solver.get_single_result(name)
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
                    except Exception:
                        pass  # Skip outputs that can't be extracted
                
                # Create DataFrames
                results_df = pd.DataFrame(results_data)
                summary_df = pd.DataFrame(summary_stats)
                
                return simulated_values, results_df, summary_df
            except Exception as e:
                print(f"Warning: Could not extract full results: {e}")
                return simulated_values, None, None
        
        return simulated_values
    
    def _iteration_callback(self, history_entry: Dict) -> None:
        """
        Callback function called after each iteration.
        
        Args:
            history_entry: Dictionary with iteration, objective, and parameters
        """
        iteration = history_entry['iteration']
        obj_value = history_entry['objective']
        params = history_entry['parameters']
        
        # Print iteration progress
        print(f"Iteration {iteration:3d}: Objective = {obj_value:.6e}", end="")
        
        # Show parameter values if not too many
        if len(params) <= 3:
            param_str = ", ".join([f"{name}={val:.3e}" for name, val in params.items()])
            print(f" | {param_str}")
        else:
            print()
    
    def _objective_function(self, param_values: np.ndarray) -> float:
        """
        Objective function wrapper for optimization.
        
        Args:
            param_values: Array of parameter values
            
        Returns:
            Objective function value
        """
        try:
            # Run simulation
            simulated_values = self._run_simulation(param_values)
            
            # Compute objective
            obj_value = self.objective_func.compute_error(simulated_values)
            
            return obj_value
        
        except Exception as e:
            # Return large value if simulation fails
            print(f"Warning: Simulation failed: {e}")
            return float(1e10)
    
    def optimize(self) -> Dict:
        """
        Run optimization.
        
        Returns:
            Dictionary with optimization results
        """
        print("Starting sv0D parameter optimization...")
        print(f"Parameters to optimize: {[p['name'] for p in self.parameters]}")
        print(f"Targets: {[t['name'] for t in self.targets]}")
        print(f"Algorithm: {self.optimization_config['algorithm']}")
        print(f"Max iterations: {self.optimization_config['max_iterations']}")
        print()
        
        # Prepare optimization inputs
        param_names = [p['name'] for p in self.parameters]
        bounds = [tuple(float(b) for b in p['bounds']) for p in self.parameters]
        
        # Get initial parameter values
        x0 = np.array([self.param_handler.get_parameter(name) for name in param_names])
        
        # Start timing
        start_time = time.time()
        
        # Run optimization with progress callback, allow graceful interruption
        interrupted = False
        result = None
        try:
            result = self.optimizer.optimize(
                objective_func=self._objective_function,
                param_names=param_names,
                bounds=bounds,
                x0=x0,
                iteration_callback=self._iteration_callback
            )
        except KeyboardInterrupt:
            interrupted = True
            print("\n\n" + "="*70)
            print("OPTIMIZATION INTERRUPTED BY USER (Ctrl+C)")
            print("="*70)
            print("Processing best results found so far...")
            print()
        
        # End timing
        end_time = time.time()
        total_time = end_time - start_time
        
        # Get results
        self.history = self.optimizer.get_history()
        self.best_value, self.best_params = self.optimizer.get_best()
        
        # For parallel differential_evolution, history tracking doesn't work properly
        # Use the result object to get actual number of function evaluations
        if result is not None:
            n_evaluations = len(self.history) if self.history else getattr(result, 'nfev', 0)
            
            # If history is empty but we have result params, use them
            if not self.history and hasattr(result, 'x') and hasattr(result, 'fun'):
                self.best_params = result.x
                self.best_value = result.fun
                if not interrupted:
                    print(f"\nNote: Detailed history not available with parallel differential_evolution")
                    print(f"Final objective value: {self.best_value:.6e}")
        else:
            # Interrupted - use best from history
            n_evaluations = len(self.history)
            if interrupted and self.best_value is not None:
                print(f"Best objective value found: {self.best_value:.6e}")
                print(f"Total evaluations completed: {n_evaluations}")
        
        # Calculate timing statistics
        avg_time_per_eval = total_time / n_evaluations if n_evaluations > 0 else 0
        
        timing_info = {
            'total_time_seconds': total_time,
            'total_time_formatted': self._format_time(total_time),
            'n_evaluations': n_evaluations,
            'avg_time_per_evaluation_seconds': avg_time_per_eval,
            'avg_time_per_evaluation_formatted': self._format_time(avg_time_per_eval)
        }
        
        # Print timing summary
        print(f"\n{'='*70}")
        print(f"TIMING SUMMARY")
        print(f"{'='*70}")
        print(f"Total optimization time: {timing_info['total_time_formatted']}")
        print(f"Total function evaluations: {n_evaluations}")
        print(f"Average time per evaluation: {timing_info['avg_time_per_evaluation_formatted']}")
        print(f"{'='*70}")
        
        # Check if optimization succeeded
        if self.best_params is None:
            if interrupted:
                error_msg = (
                    "Cannot process results: optimization was interrupted too early.\n"
                    "When using parallel differential_evolution, history tracking is not available.\n"
                    "Try interrupting later, or set 'parallel: false' in tuning.yaml for full history tracking."
                )
            else:
                error_msg = "Optimization failed: no best parameters found"
            raise RuntimeError(error_msg)
        
        # Run final simulation with best parameters and get full results
        print("\nRunning final simulation with best parameters...")
        should_extract_full = self.output_config.get('save_plots', True)
        
        if should_extract_full:
            final_simulated_values, optimized_results_df, optimized_summary_df = \
                self._run_simulation(self.best_params, return_full_results=True)
        else:
            final_simulated_values = self._run_simulation(self.best_params)
            optimized_results_df = None
            optimized_summary_df = None
        
        # Save results
        if self.output_config.get('save_final_config', True):
            # Generate output filename based on original model config filename
            original_filename = os.path.basename(self.config['model']['config_file'])
            name_without_ext, ext = os.path.splitext(original_filename)
            optimized_filename = f"{name_without_ext}_optimized{ext}"
            
            final_config = self.param_handler.get_config()
            self.result_handler.save_final_config(
                final_config,
                filename=optimized_filename
            )
        
        # Generate report
        self.result_handler.create_summary_report(
            history=self.history,
            param_names=param_names,
            best_value=self.best_value,
            best_params=self.best_params,
            targets=self.targets,
            simulated_values=final_simulated_values,
            optimized_results_df=optimized_results_df,
            optimized_summary_df=optimized_summary_df,
            timing_info=timing_info
        )
        
        return {
            'success': result.success if result else False,
            'message': result.message if result else ('Interrupted by user' if interrupted else 'Unknown error'),
            'best_value': self.best_value,
            'best_params': dict(zip(param_names, self.best_params)),
            'history': self.history,
            'result': result,
            'interrupted': interrupted
        }
    
    def evaluate(self, param_values: Optional[Dict[str, float]] = None) -> Dict:
        """
        Evaluate objective function without optimization.
        
        Args:
            param_values: Dictionary of parameter values (optional, uses current if None)
            
        Returns:
            Dictionary with evaluation results
        """
        if param_values:
            # Update parameters
            for name, value in param_values.items():
                self.param_handler.set_parameter(name, value)
        
        # Get current parameter values
        param_names = [p['name'] for p in self.parameters]
        current_values = np.array([
            self.param_handler.get_parameter(name) for name in param_names
        ])
        
        # Run simulation
        simulated_values = self._run_simulation(current_values)
        
        # Compute objective
        obj_value = self.objective_func.compute_error(simulated_values)
        
        return {
            'objective_value': obj_value,
            'simulated_values': simulated_values,
            'parameters': dict(zip(param_names, current_values))
        }


def run_optimization(config_file: str) -> Dict:
    """
    Convenience function to run optimization from config file.
    
    Args:
        config_file: Path to YAML configuration file
        
    Returns:
        Optimization results dictionary
    """
    tuner = SV0DTuner(config_file)
    return tuner.optimize()
