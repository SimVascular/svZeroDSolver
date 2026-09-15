"""
Sensitivity Analysis for sv0D Tuning Framework

Performs correlation-based sensitivity screening to understand how parameter
variations affect quantities of interest.

This implementation uses Sobol low-discrepancy sampling for parameter space
exploration and reports screening metrics rather than true Sobol indices.
The method provides screening-level sensitivity information suitable for
identifying the most influential parameters.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Any, Callable
from scipy.stats import qmc
import json
import os
from pathlib import Path
import yaml

from .parameter_handler import ParameterHandler
from .output_extractor import OutputExtractor
from .simulation import run_simulation
from .expression_handler import Expression


class SensitivityAnalyzer:
    """
    Performs correlation-based sensitivity screening.
    """

    ANALYSIS_LABEL = "correlation-based sensitivity screening"
    FIRST_ORDER_LABEL = "First-order screening score (squared correlation)"
    TOTAL_ORDER_LABEL = (
        "Total-order screening score (binned conditional-variance heuristic)"
    )
    INDICES_FILENAME = "screening_indices.json"
    PER_QOI_FIGURE_PREFIX = "screening"
    HEATMAP_FIRST_FILENAME = "screening_heatmap_first_order.png"
    HEATMAP_TOTAL_FILENAME = "screening_heatmap_total_order.png"
    
    def __init__(self, config_file: str):
        """
        Initialize sensitivity analyzer with configuration file.
        
        Args:
            config_file: Path to sensitivity configuration YAML file
        """
        # Load configuration directly (don't use ConfigHandler which expects optimization format)
        self.config_file = config_file
        self.config = self._load_config()
        self._validate_config()
        
        # Get model config file path
        config_dir = os.path.dirname(os.path.abspath(self.config_file))
        model_config = self.config['model']['config_file']
        if not os.path.isabs(model_config):
            model_config = os.path.join(config_dir, model_config)
        
        # Initialize parameter handler
        self.param_handler = ParameterHandler(model_config)
        
        # Get configuration sections
        self.parameters = self.config.get('parameters', [])
        self.quantities_of_interest = self.config.get('quantities_of_interest', [])
        self.sensitivity_config = self.config.get('sensitivity', {})
        self.output_config = self.config.get('output', {})

        for param in self.parameters:
            try:
                self.param_handler.get_parameter(param["name"])
            except ValueError as exc:
                raise ValueError(
                    f"Sensitivity parameter '{param['name']}' not found in model configuration"
                ) from exc
        
        # Sensitivity analysis settings
        self.n_samples = self.sensitivity_config.get('n_samples', 512)

        # Replace expression string with Expression object for each QoI (all scalars)
        for qoi in self.quantities_of_interest:
            expr_str = qoi["expression"]
            qoi["expression"] = Expression(expr_str, "scalar")

        # Results storage
        self.results = {}
        self.sample_data = []
    
    def _load_config(self) -> Dict:
        """Load YAML configuration file."""
        with open(self.config_file, 'r') as f:
            return yaml.safe_load(f)
    
    def _validate_config(self):
        """Validate configuration structure."""
        required_sections = ['model', 'parameters', 'quantities_of_interest']
        for section in required_sections:
            if section not in self.config:
                raise ValueError(f"Missing required section '{section}' in sensitivity configuration")
        
        # Validate model section
        if 'config_file' not in self.config['model']:
            raise ValueError("model.config_file is required")
        
        # Validate parameters
        if not self.config['parameters']:
            raise ValueError("At least one parameter must be specified")
        
        for i, param in enumerate(self.config['parameters']):
            if 'name' not in param:
                raise ValueError(f"Parameter {i} missing 'name'")
            if 'bounds' not in param:
                raise ValueError(f"Parameter '{param['name']}' missing 'bounds'")
            if len(param['bounds']) != 2:
                raise ValueError(f"Parameter '{param['name']}' bounds must have 2 values [min, max]")
        
        # Validate quantities of interest
        if not self.config['quantities_of_interest']:
            raise ValueError("At least one quantity of interest must be specified")
        
        for i, qoi in enumerate(self.config['quantities_of_interest']):
            if 'name' not in qoi:
                raise ValueError(f"Quantity of interest {i} missing 'name'")
            if 'expression' not in qoi:
                raise ValueError(
                    f"Quantity of interest '{qoi.get('name', i)}' missing 'expression'"
                )
        
    def _get_simulated_quantities_of_interest(
        self, param_values: np.ndarray, suppress_warnings: bool = False
    ) -> Dict[str, float]:
        """
        Run sv0D simulation and return quantity-of-interest values.
        Each QoI has name (free label) and expression (e.g. np.max(pressure:AV:AR_SYS)).
        """
        qoi_values = {qoi["name"]: np.nan for qoi in self.quantities_of_interest}

        try:
            _, extractor = run_simulation(
                self.param_handler, self.parameters, param_values
            )
            times = extractor.get_times()
            available_outputs = extractor.get_all_output_names()

            # Collect outputs needed by any QoI expression
            unique_outputs = {}
            for qoi in self.quantities_of_interest:
                expr = qoi["expression"]
                for out_name in expr.output_names(available_outputs):
                    if out_name not in unique_outputs:
                        try:
                            ts = extractor.extract(out_name, "time_series")
                            unique_outputs[out_name] = {
                                "time_series": ts,
                                "times": times,
                            }
                        except Exception:
                            pass

            for qoi in self.quantities_of_interest:
                name = qoi["name"]
                try:
                    qoi_values[name] = float(
                        qoi["expression"].evaluate(
                            unique_outputs, available_outputs
                        )
                    )
                except Exception as e:
                    if not suppress_warnings:
                        print(f"Warning: Could not evaluate QoI '{name}': {e}")
                    qoi_values[name] = np.nan

        except Exception as e:
            if not suppress_warnings:
                print(f"Warning: Simulation failed: {e}")

        return qoi_values
    
    def _create_qoi_function(self, qoi_key: str) -> Callable:
        """
        Create a function that evaluates a specific QoI for given parameters.
        
        Args:
            qoi_key: Key for the quantity of interest
            
        Returns:
            Function that takes parameter array and returns QoI value
        """
        def qoi_func(params):
            """Evaluate QoI for given parameters."""
            try:
                # Handle both single parameter set and multiple parameter sets
                if params.ndim == 1:
                    # Single parameter set
                    qoi_values = self._get_simulated_quantities_of_interest(params)
                    return qoi_values[qoi_key]
                else:
                    # Multiple parameter sets
                    results = []
                    for i in range(params.shape[0]):
                        qoi_values = self._get_simulated_quantities_of_interest(params[i, :])
                        results.append(qoi_values[qoi_key])
                    return np.array(results)
            except Exception as e:
                print(f"Error in simulation: {e}")
                if params.ndim == 1:
                    return np.nan
                else:
                    return np.full(params.shape[0], np.nan)
        
        return qoi_func
    
    def run(self) -> Dict:
        """
        Run correlation-based sensitivity screening.
        
        Returns:
            Dictionary with sensitivity analysis results
        """
        print("="*70)
        print("SENSITIVITY ANALYSIS")
        print("="*70)
        print(f"Parameters: {[p['name'] for p in self.parameters]}")
        qoi_names = [q["name"] for q in self.quantities_of_interest]
        print(f"Quantities of Interest: {qoi_names}")
        print(f"Number of samples: {self.n_samples}")
        print()

        # Reset per-run state so reusing the analyzer does not mix old and new samples.
        self.results = {}
        self.sample_data = []
        
        # Get parameter bounds
        param_names = [p['name'] for p in self.parameters]
        bounds = [p['bounds'] for p in self.parameters]
        n_params = len(param_names)
        
        # Generate Sobol low-discrepancy samples for screening.
        print("Generating quasi-random screening samples...")
        
        # Extract and validate bounds
        lower_bounds = []
        upper_bounds = []
        for i, (param, bound) in enumerate(zip(self.parameters, bounds)):
            lower = float(bound[0])
            upper = float(bound[1])
            if lower >= upper:
                raise ValueError(
                    f"Invalid bounds for parameter '{param['name']}': "
                    f"lower ({lower}) must be < upper ({upper})"
                )
            lower_bounds.append(lower)
            upper_bounds.append(upper)
            print(f"  {param['name']}: [{lower:.3e}, {upper:.3e}]")
        
        lower_bounds = np.array(lower_bounds)
        upper_bounds = np.array(upper_bounds)
        
        sampler = qmc.Sobol(d=n_params, scramble=True)
        # Generate samples in [0, 1]^d
        samples_unit = sampler.random(self.n_samples)
        
        # Scale to parameter bounds
        samples = qmc.scale(samples_unit, lower_bounds, upper_bounds)
        
        print(f"Generated {samples.shape[0]} samples in {n_params}-dimensional space")
        print()
        
        # Evaluate all QoIs for all samples
        print("Evaluating simulations...")
        print("-"*70)
        
        all_qoi_values = {qoi["name"]: [] for qoi in self.quantities_of_interest}
        
        for i, param_values in enumerate(samples):
            if (i + 1) % max(1, self.n_samples // 20) == 0:
                print(f"  Progress: {i+1}/{self.n_samples} ({100*(i+1)/self.n_samples:.0f}%)")
            
            qoi_values = self._get_simulated_quantities_of_interest(param_values)
            
            # Store sample data
            sample_entry = {
                'sample_id': i,
                **{param_names[j]: param_values[j] for j in range(n_params)},
                **qoi_values
            }
            self.sample_data.append(sample_entry)
            
            # Store QoI values
            for qoi_key, qoi_value in qoi_values.items():
                all_qoi_values[qoi_key].append(qoi_value)
        
        print(f"  Completed {self.n_samples} simulations")
        print()
        
        # Compute screening metrics for each QoI
        print("Computing screening metrics...")
        print("-"*70)
        
        for qoi_key in all_qoi_values.keys():
            print(f"  Analyzing {qoi_key}...")
            
            qoi_values_array = np.array(all_qoi_values[qoi_key])
            
            # Check for NaN values
            if np.any(np.isnan(qoi_values_array)):
                n_nan = np.sum(np.isnan(qoi_values_array))
                print(f"    Warning: {n_nan}/{len(qoi_values_array)} simulations failed")
                # Remove NaN values for analysis
                valid_idx = ~np.isnan(qoi_values_array)
                qoi_values_array = qoi_values_array[valid_idx]
                samples_valid = samples[valid_idx, :]
            else:
                samples_valid = samples
            
            # Compute screening metrics from the sampled outputs.
            try:
                total_variance = np.var(qoi_values_array)
                
                if total_variance < 1e-12:
                    print(f"    Warning: Very low variance ({total_variance:.3e}), QoI may be insensitive to parameters")
                    # All parameters have zero sensitivity
                    first_order = {param_names[i]: 0.0 for i in range(n_params)}
                    total_order = first_order.copy()
                else:
                    first_order = {}
                    total_order = {}
                    
                    # Estimate first-order and total-order screening scores.
                    for i, param_name in enumerate(param_names):
                        param_values = samples_valid[:, i]
                        
                        # First-order score: squared correlation coefficient
                        correlation = np.corrcoef(param_values, qoi_values_array)[0, 1]
                        first_order[param_name] = max(0.0, min(1.0, correlation**2))
                        
                        # Total-order score: use conditional variance estimation
                        # Group samples by parameter value bins
                        n_bins = min(10, len(param_values) // 10)
                        if n_bins >= 3:
                            bins = np.percentile(param_values, np.linspace(0, 100, n_bins+1))
                            bin_indices = np.digitize(param_values, bins[1:-1])
                            
                            # Variance within bins (conditional variance)
                            var_within = 0.0
                            for bin_idx in range(n_bins):
                                mask = bin_indices == bin_idx
                                if np.sum(mask) > 1:
                                    var_within += np.var(qoi_values_array[mask]) * np.sum(mask)
                            var_within /= len(qoi_values_array)
                            
                            # Total-order screening score: 1 - conditional variance fraction
                            total_order[param_name] = max(0.0, min(1.0, 1.0 - var_within / total_variance))
                        else:
                            # Not enough samples for binning, use first-order as approximation
                            total_order[param_name] = first_order[param_name]
                
                self.results[qoi_key] = {
                    'analysis_type': 'correlation_screening',
                    'sampler': 'sobol_sequence',
                    'first_order_metric': 'squared_pearson_correlation',
                    'total_order_metric': 'binned_conditional_variance_screening',
                    'first_order': first_order,
                    'total_order': total_order,
                    'mean': float(np.mean(qoi_values_array)),
                    'std': float(np.std(qoi_values_array)),
                    'min': float(np.min(qoi_values_array)),
                    'max': float(np.max(qoi_values_array))
                }
                
            except Exception as e:
                print(f"    Error computing indices: {e}")
                import traceback
                traceback.print_exc()
                # Fallback to zeros
                self.results[qoi_key] = {
                    'analysis_type': 'correlation_screening',
                    'sampler': 'sobol_sequence',
                    'first_order_metric': 'squared_pearson_correlation',
                    'total_order_metric': 'binned_conditional_variance_screening',
                    'first_order': {param_names[i]: 0.0 for i in range(n_params)},
                    'total_order': {param_names[i]: 0.0 for i in range(n_params)},
                    'mean': float(np.mean(qoi_values_array)),
                    'std': float(np.std(qoi_values_array)),
                    'min': float(np.min(qoi_values_array)),
                    'max': float(np.max(qoi_values_array))
                }
        
        print()
        print("✓ Sensitivity analysis complete")
        print()
        
        return self.results
    
    def save_results(self, output_dir: Optional[str] = None):
        """
        Save sensitivity analysis results to files.
        
        Args:
            output_dir: Output directory (uses config if not provided)
        """
        if output_dir is None:
            output_dir = self.output_config.get('directory', 'sensitivity_results')
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        print(f"Saving results to {output_dir}/")
        print("-"*70)
        
        # Save screening metrics as JSON
        indices_file = output_path / self.INDICES_FILENAME
        with open(indices_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"  ✓ Saved screening metrics: {indices_file.name}")
        
        # Save sample data as CSV
        if self.sample_data:
            samples_df = pd.DataFrame(self.sample_data)
            samples_file = output_path / 'sample_data.csv'
            samples_df.to_csv(samples_file, index=False)
            print(f"  ✓ Saved sample data: {samples_file.name}")
        
        # Create summary report
        self._create_summary_report(output_path)
        print(f"  ✓ Saved summary report: summary.txt")
        
        # Create visualizations
        try:
            self._create_visualizations(output_path)
            print(f"  ✓ Saved bar charts and scatter plots")
        except Exception as e:
            print(f"  Warning: Could not create visualizations: {e}")
        
        print()
    
    def _create_summary_report(self, output_path: Path):
        """Create a text summary report."""
        report_file = output_path / 'summary.txt'
        
        with open(report_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("SENSITIVITY ANALYSIS SUMMARY\n")
            f.write("="*70 + "\n\n")
            f.write(f"Method: {self.ANALYSIS_LABEL}\n")
            f.write("Sampling: Sobol low-discrepancy sequence\n\n")
            
            f.write(f"Parameters analyzed: {[p['name'] for p in self.parameters]}\n")
            f.write(f"Number of samples: {self.n_samples}\n\n")
            
            for qoi_key, results in self.results.items():
                f.write("-"*70 + "\n")
                f.write(f"Quantity of Interest: {qoi_key}\n")
                f.write("-"*70 + "\n")
                f.write(f"Mean: {results['mean']:.6e}\n")
                f.write(f"Std:  {results['std']:.6e}\n")
                f.write(f"Min:  {results['min']:.6e}\n")
                f.write(f"Max:  {results['max']:.6e}\n\n")
                
                f.write(f"{self.FIRST_ORDER_LABEL}:\n")
                for param, value in results['first_order'].items():
                    f.write(f"  {param:<30} {value:>10.4f}\n")
                f.write("\n")
                
                f.write(f"{self.TOTAL_ORDER_LABEL}:\n")
                for param, value in results['total_order'].items():
                    f.write(f"  {param:<30} {value:>10.4f}\n")
                f.write("\n\n")
    
    def _create_visualizations(self, output_path: Path):
        """Create visualization plots."""
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib not available, skipping visualizations")
            return
        
        param_names = [p['name'] for p in self.parameters]
        
        # Create heatmap of all screening scores
        self._create_heatmap(output_path, param_names)
        
        for qoi_key, results in self.results.items():
            # Create bar plot of screening scores
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
            
            # First-order scores
            first_order_values = [results['first_order'][p] for p in param_names]
            ax1.bar(range(len(param_names)), first_order_values)
            ax1.set_xticks(range(len(param_names)))
            ax1.set_xticklabels(param_names, rotation=45, ha='right')
            ax1.set_ylabel('First-order screening score')
            ax1.set_title('Main Effects')
            ax1.grid(axis='y', alpha=0.3)
            ax1.set_ylim([0, 1])
            
            # Total-order scores
            total_order_values = [results['total_order'][p] for p in param_names]
            ax2.bar(range(len(param_names)), total_order_values)
            ax2.set_xticks(range(len(param_names)))
            ax2.set_xticklabels(param_names, rotation=45, ha='right')
            ax2.set_ylabel('Total-order screening score')
            ax2.set_title('Total Effects (Main + Interactions)')
            ax2.grid(axis='y', alpha=0.3)
            ax2.set_ylim([0, 1])
            
            plt.suptitle(f'Sensitivity Screening: {qoi_key}')
            plt.tight_layout()
            
            # Save figure
            safe_filename = qoi_key.replace(':', '_').replace('/', '_')
            fig_file = output_path / f'{self.PER_QOI_FIGURE_PREFIX}_{safe_filename}.png'
            plt.savefig(fig_file, dpi=150, bbox_inches='tight')
            plt.close()
        
        # Create scatter plots of QoI vs parameters
        if self.sample_data:
            samples_df = pd.DataFrame(self.sample_data)
            
            for qoi_key in self.results.keys():
                if qoi_key not in samples_df.columns:
                    continue
                
                n_params = len(param_names)
                fig, axes = plt.subplots(1, n_params, figsize=(5*n_params, 4))
                if n_params == 1:
                    axes = [axes]
                
                for i, param_name in enumerate(param_names):
                    axes[i].scatter(samples_df[param_name], samples_df[qoi_key], 
                                   alpha=0.5, s=10)
                    axes[i].set_xlabel(param_name)
                    axes[i].set_ylabel(qoi_key)
                    axes[i].grid(alpha=0.3)
                
                plt.suptitle(f'Parameter Effects on {qoi_key}')
                plt.tight_layout()
                
                safe_filename = qoi_key.replace(':', '_').replace('/', '_')
                fig_file = output_path / f'scatter_{safe_filename}.png'
                plt.savefig(fig_file, dpi=150, bbox_inches='tight')
                plt.close()
    
    def _create_heatmap(self, output_path: Path, param_names: List[str]):
        """Create heatmap visualization of all screening scores.
        
        Rows = parameters, Columns = quantities of interest.
        """
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            return
        
        qoi_keys = list(self.results.keys())
        
        # Create matrices: rows = parameters, columns = quantities of interest
        first_order_matrix = np.zeros((len(param_names), len(qoi_keys)))
        total_order_matrix = np.zeros((len(param_names), len(qoi_keys)))
        
        for j, param_name in enumerate(param_names):
            for i, qoi_key in enumerate(qoi_keys):
                first_order_matrix[j, i] = self.results[qoi_key]['first_order'][param_name]
                total_order_matrix[j, i] = self.results[qoi_key]['total_order'][param_name]
        
        figsize = (4 + len(qoi_keys)*1.2, len(param_names)*0.8 + 2)
        
        # Figure 1: First-order heatmap
        fig1, ax1 = plt.subplots(1, 1, figsize=figsize)
        im1 = ax1.imshow(first_order_matrix, aspect='auto', cmap='YlOrRd', vmin=0, vmax=1)
        ax1.set_xticks(range(len(qoi_keys)))
        ax1.set_yticks(range(len(param_names)))
        ax1.set_xticklabels(qoi_keys, rotation=45, ha='right')
        ax1.set_yticklabels(param_names, fontsize=9)
        ax1.set_title(self.FIRST_ORDER_LABEL, fontweight='bold')
        ax1.set_xlabel('Quantities of Interest')
        ax1.set_ylabel('Parameters')
        for j in range(len(param_names)):
            for i in range(len(qoi_keys)):
                value = first_order_matrix[j, i]
                color = 'white' if value > 0.5 else 'black'
                ax1.text(i, j, f'{value:.3f}', ha='center', va='center', 
                        color=color, fontsize=8, fontweight='bold')
        cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
        cbar1.set_label('Sensitivity Index', rotation=270, labelpad=15)
        plt.tight_layout()
        fig1.savefig(output_path / self.HEATMAP_FIRST_FILENAME, dpi=200, bbox_inches='tight')
        plt.close(fig1)
        
        # Figure 2: Total-order heatmap
        fig2, ax2 = plt.subplots(1, 1, figsize=figsize)
        im2 = ax2.imshow(total_order_matrix, aspect='auto', cmap='YlOrRd', vmin=0, vmax=1)
        ax2.set_xticks(range(len(qoi_keys)))
        ax2.set_yticks(range(len(param_names)))
        ax2.set_xticklabels(qoi_keys, rotation=45, ha='right')
        ax2.set_yticklabels(param_names, fontsize=9)
        ax2.set_title(self.TOTAL_ORDER_LABEL, fontweight='bold')
        ax2.set_xlabel('Quantities of Interest')
        ax2.set_ylabel('Parameters')
        for j in range(len(param_names)):
            for i in range(len(qoi_keys)):
                value = total_order_matrix[j, i]
                color = 'white' if value > 0.5 else 'black'
                ax2.text(i, j, f'{value:.3f}', ha='center', va='center', 
                        color=color, fontsize=8, fontweight='bold')
        cbar2 = plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
        cbar2.set_label('Sensitivity Index', rotation=270, labelpad=15)
        plt.tight_layout()
        fig2.savefig(output_path / self.HEATMAP_TOTAL_FILENAME, dpi=200, bbox_inches='tight')
        plt.close(fig2)
        
        print(
            f"  ✓ Created heatmaps: {self.HEATMAP_FIRST_FILENAME}, {self.HEATMAP_TOTAL_FILENAME}"
        )


def run_sensitivity_analysis(config_file: str) -> Dict:
    """
    Convenience function to run sensitivity analysis from config file.
    
    Args:
        config_file: Path to sensitivity configuration YAML file
        
    Returns:
        Sensitivity analysis results dictionary
    """
    analyzer = SensitivityAnalyzer(config_file)
    results = analyzer.run()
    analyzer.save_results()
    return results
