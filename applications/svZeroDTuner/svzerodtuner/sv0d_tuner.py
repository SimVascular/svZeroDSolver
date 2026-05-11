"""
Main sv0D Tuning Framework

Orchestrates parameter optimization to match target output values.
"""

import os
import json
import time
import numpy as np
import pandas as pd
import multiprocessing as mp
from functools import partial
from typing import Dict, List, Optional, Callable

from .parameter_handler import ParameterHandler
from .output_extractor import OutputExtractor
from .simulation import run_simulation
from .objective import create_objective, ObjectiveFunction
from .optimizer import OptimizerWrapper
from .config_handler import ConfigHandler
from .result_handler import ResultHandler
from .scaling import get_scaling, to_physical_array, to_opt_array

from .expression_handler import Expression


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
        self.objective_func = create_objective(targets=self.targets, **self.objective_config)

        # Replace expression string with Expression object for each target
        for target in self.targets:
            expr_str = target.get("expression")
            if expr_str:
                target["expression"] = Expression(
                    expr_str, target.get("type", "time_series")
                )

        # Initialize scaling objects for each parameter (optimizer space <-> physical space)
        self._scalings = [
            get_scaling(p.get("scaling", "identity"), tuple(p["bounds"]) if "bounds" in p else None)
            for p in self.parameters
        ]

        # Validate all parameter names exist in the model at init time
        for p in self.parameters:
            try:
                self.param_handler.get_parameter(p["name"])
            except ValueError as e:
                raise ValueError(f"Invalid parameter in tuning config: {e}") from e

        # State
        self.solver = None
        self.extractor = None
        self.history = []
        self.best_value = None
        self.best_params = None

    def __getstate__(self):
        """
        Make tuner pickle-safe for multiprocessing objective dispatch.
        `pysvzerod.Solver`/extractor instances are runtime-only and not pickleable.
        """
        state = self.__dict__.copy()
        state["solver"] = None
        state["extractor"] = None
        return state
    
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
    
    def _get_simulated_values(
        self, param_values: np.ndarray, return_full_results: bool = False
    ):
        """
        Run sv0D simulation and return simulated target values.

        Args:
            param_values: Array of parameter values
            return_full_results: If True, return tuple of (simulated_values, results_df, summary_df)

        Returns:
            If return_full_results=False: Dictionary of simulated output values
            If return_full_results=True: Tuple of (simulated_values dict, results_df, summary_df)
        """
        self.solver, self.extractor = run_simulation(
            self.param_handler, self.parameters, param_values
        )

        # Extract all target outputs; store only simulated_values[name] per target
        simulated_values = {}
        times = self.extractor.get_times()
        available_outputs = self.extractor.get_all_output_names()

        # Collect output names we need from each target's expression
        unique_outputs = {}
        for target in self.targets:
            expr = target.get("expression")
            if not expr or not isinstance(expr, Expression):
                raise ValueError(f"Target '{target['name']}' must have an Expression object")
            for out_name in expr.output_names(available_outputs):
                if out_name not in unique_outputs:
                    try:
                        ts = self.extractor.extract(out_name, "time_series")
                        unique_outputs[out_name] = {
                            "time_series": ts,
                            "times": times,
                        }
                    except Exception:
                        raise ValueError(f"Failed to extract output '{out_name}' for target '{target['name']}'")

        # Compute simulated value per target; store under target name
        for target in self.targets:
            name = target["name"]
            expr = target.get("expression")
            try:
                result = expr.evaluate(unique_outputs, available_outputs)
                if expr.kind == "time_series":
                    arr, t = result
                    simulated_values[name] = arr
                    simulated_values[f"{name}_times"] = t
                else:
                    simulated_values[name] = float(result)
            except Exception as e:
                raise RuntimeError(
                    f"Failed to evaluate target '{name}': {e}"
                ) from e
        
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
            simulated_values = self._get_simulated_values(param_values)
            
            # Compute objective
            obj_value = self.objective_func.compute(simulated_values)
            
            return obj_value
        
        except Exception as e:
            # Re-raise config errors; only penalize numerical solver failures
            if isinstance(e, (ValueError, KeyError, AttributeError)):
                raise
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
        print(f"Optimizer configuration: {self.optimization_config}")
        if self.optimization_config['algorithm'] == 'differential_evolution':
            n_workers = self.optimization_config.get('workers', 1)
            if n_workers == -1:
                n_workers = mp.cpu_count()
            print(f"Running in parallel with {n_workers} workers")
        print()
        
        # Prepare optimization inputs
        param_names = [p['name'] for p in self.parameters]
        bounds = [tuple(float(b) for b in p['bounds']) for p in self.parameters]
        
        # Get initial parameter values
        print(f"Using initial parameter values from JSON file:")
        for name in param_names:
            print(f"\t{name}: {self.param_handler.get_parameter(name)}")
        x0 = np.array([self.param_handler.get_parameter(name) for name in param_names])

        # Create scaling functions for optimizer space <-> physical space
        to_opt = partial(to_opt_array, scalings=self._scalings)
        to_phys = partial(to_physical_array, scalings=self._scalings)
        # Start timing
        start_time = time.time()
        
        # Run optimization, allow graceful interruption
        interrupted = False
        result = None
        try:
            result = self.optimizer.optimize(
                objective_func=self._objective_function,
                param_names=param_names,
                bounds=bounds,
                x0=x0,
                parameters=self.parameters,
                param_scaling_to_opt_space=to_opt,
                param_scaling_to_phys_space=to_phys,
            )
        except KeyboardInterrupt:
            interrupted = True
            print("\n\n" + "="*70)
            print("OPTIMIZATION INTERRUPTED BY USER (Ctrl+C)")
            print("="*70)
            print("Processing best results found so far...")
            print()
        
        # Print termination reason when optimization completed normally
        if result is not None and not interrupted:
            print("\n" + "="*70)
            print("OPTIMIZATION TERMINATION")
            print("="*70)
            success = getattr(result, 'success', None)
            message = getattr(result, 'message', None)
            status = getattr(result, 'status', None)
            if message:
                print(f"Reason: {message}")
            if success is not None:
                print(f"Success: {success}")
            if status is not None:
                print(f"Status code: {status}")
            nfev = getattr(result, 'nfev', None)
            if nfev is not None:
                print(f"Function evaluations: {nfev}")
            nit = getattr(result, 'nit', None)
            if nit is not None:
                print(f"Algorithm iterations: {nit}")
            print("="*70)
        
        # End timing
        end_time = time.time()
        total_time = end_time - start_time
        
        # Get results
        self.history = self.optimizer.get_history()
        self.best_value, self.best_params = self.optimizer.get_best()

        if result is not None:
            n_evaluations = getattr(result, 'nfev', len(self.history) if self.history else 0)
            n_iterations = getattr(result, 'nit', len(self.history) if self.history else 0)
        else:
            # Interrupted - use best from history; nfev unknown
            n_evaluations = None  # not available when interrupted
            n_iterations = len(self.history)
            if interrupted and self.best_value is not None:
                print(f"Best objective value found: {self.best_value:.6e}")
                print(f"Total iterations completed: {n_iterations}")

        # Calculate timing statistics
        avg_time_per_eval = total_time / n_evaluations if n_evaluations and n_evaluations > 0 else 0
        timing_info = {
            'total_time_seconds': total_time,
            'total_time_formatted': self._format_time(total_time),
            'n_evaluations': n_evaluations,
            'avg_time_per_evaluation_seconds': avg_time_per_eval,
            'avg_time_per_evaluation_formatted': self._format_time(avg_time_per_eval)
        }
        avg_time_per_iter = total_time / n_iterations if n_iterations > 0 else 0
        timing_info.update({
            'n_iterations': n_iterations,
            'avg_time_per_iteration_seconds': avg_time_per_iter,
            'avg_time_per_iteration_formatted': self._format_time(avg_time_per_iter)
        })

        # Print timing summary
        print(f"\n{'='*70}")
        print(f"TIMING SUMMARY")
        print(f"{'='*70}")
        print(f"Total optimization time: {timing_info['total_time_formatted']}")
        print(f"Total iterations: {n_iterations}")
        print(f"Total function evaluations: {n_evaluations if n_evaluations is not None else 'not available (interrupted)'}")
        print(f"Average time per iteration: {timing_info['avg_time_per_iteration_formatted']}")
        if n_evaluations is not None:
            print(f"Average time per evaluation: {timing_info['avg_time_per_evaluation_formatted']}")
        else:
            print(f"Average time per evaluation: not available (interrupted)")
        print(f"{'='*70}")
        
        # Check if optimization succeeded
        if self.best_params is None:
            if interrupted:
                error_msg = (
                    "Cannot process results: optimization was interrupted before any results were obtained.\n"
                    "Try interrupting later to allow more progress."
                )
            else:
                error_msg = "Optimization failed: no best parameters found"
            raise RuntimeError(error_msg)
        
        # Run final simulation with best parameters and get full results
        print("\nRunning final simulation with best parameters...")
        should_extract_full = self.output_config.get('save_plots', True)
        
        if should_extract_full:
            final_simulated_values, optimized_results_df, optimized_summary_df = \
                self._get_simulated_values(self.best_params, return_full_results=True)
        else:
            final_simulated_values = self._get_simulated_values(self.best_params)
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
            parameters=self.parameters,
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
        simulated_values = self._get_simulated_values(current_values)
        
        # Compute objective
        obj_value = self.objective_func.compute(simulated_values)
        
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
