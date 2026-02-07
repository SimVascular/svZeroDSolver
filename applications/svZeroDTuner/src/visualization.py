"""
Visualization functions for sv0D Tuning Framework

Plots optimization history including objective value convergence and parameter evolution.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict, Optional


def plot_objective_history(
    history: List[Dict],
    output_file: Optional[str] = None,
    show: bool = False
):
    """
    Plot objective value vs function evaluation.
    
    Args:
        history: Optimization history list
        output_file: Path to save plot (optional)
        show: Whether to display plot
    """
    if not history:
        return
    
    evaluations = [h['evaluation'] for h in history]
    objectives = [h['objective'] for h in history]
    
    plt.figure(figsize=(10, 6))
    plt.plot(evaluations, objectives, 'ko-', linewidth=2, markersize=6, label='Objective Value')
    plt.xlabel('Function evaluation', fontsize=12)
    plt.ylabel('Objective Value', fontsize=12)
    plt.title('Optimization Convergence', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    
    # Add best value annotation
    best_idx = np.argmin(objectives)
    best_eval = evaluations[best_idx]
    best_obj = objectives[best_idx]
    plt.plot(best_eval, best_obj, 'ro', markersize=10, label=f'Best: {best_obj:.6e}')
    plt.annotate(
        f'Best: {best_obj:.6e}',
        xy=(best_eval, best_obj),
        xytext=(10, 10),
        textcoords='offset points',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
    )
    
    # Set y-axis limits: slightly below 0 to initial objective value (with small margin)
    initial_obj = objectives[0]
    y_min = -initial_obj * 0.02  # 2% below zero
    y_max = initial_obj * 1.05  # Add 5% margin at top
    plt.ylim(y_min, y_max)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved objective history plot to {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_parameter_evolution(
    history: List[Dict],
    param_names: List[str],
    output_file: Optional[str] = None,
    show: bool = False
):
    """
    Plot parameter evolution over function evaluations.
    
    Args:
        history: Optimization history list
        param_names: List of parameter names
        output_file: Path to save plot (optional)
        show: Whether to display plot
    """
    if not history or not param_names:
        return
    
    n_params = len(param_names)
    n_cols = min(3, n_params)
    n_rows = (n_params + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_params == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = axes if isinstance(axes, np.ndarray) else [axes]
    else:
        axes = axes.flatten()
    
    evaluations = [h['evaluation'] for h in history]
    
    for idx, param_name in enumerate(param_names):
        ax = axes[idx]
        param_values = [h['parameters'][param_name] for h in history]
        
        ax.plot(evaluations, param_values, 'ko-', linewidth=2, markersize=6)
        ax.set_xlabel('Function evaluation', fontsize=10)
        ax.set_ylabel(param_name, fontsize=10)
        ax.set_title(f'Parameter: {param_name}', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Add initial and final values
        ax.axhline(y=param_values[0], color='g', linestyle='--', alpha=0.5, label='Initial')
        ax.axhline(y=param_values[-1], color='r', linestyle='--', alpha=0.5, label='Final')
        ax.legend(fontsize=9)
    
    # Hide unused subplots
    for idx in range(n_params, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved parameter evolution plot to {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_target_comparison(
    targets: List[Dict],
    simulated_values: Dict,
    output_file: Optional[str] = None,
    show: bool = False
):
    """
    Plot target vs simulated values comparison.
    
    Args:
        targets: List of target specifications
        simulated_values: Dictionary of simulated values
        output_file: Path to save plot (optional)
        show: Whether to display plot
    """
    n_targets = len(targets)
    if n_targets == 0:
        return
    
    n_cols = min(3, n_targets)
    n_rows = (n_targets + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_targets == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = axes if isinstance(axes, np.ndarray) else [axes]
    else:
        axes = axes.flatten()
    
    for idx, target in enumerate(targets):
        ax = axes[idx]
        name = target['name']
        target_type = target.get('type', 'time_series')
        
        # Construct the key used to store the simulated value
        # For scalar targets, the key includes the type (e.g., "pressure:LV:AV_max")
        if target_type in ['min', 'max', 'mean']:
            sim_key = f"{name}_{target_type}"
        else:
            sim_key = name
        
        if sim_key not in simulated_values:
            ax.text(0.5, 0.5, f'No data for {name}', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Target: {name}', fontsize=11, fontweight='bold')
            continue
        
        sim_value = simulated_values[sim_key]
        
        if target_type == 'time_series':
            # Plot time series comparison
            if 'target_times' in target and 'target_values' in target:
                target_times = target['target_times']
                target_values = target['target_values']
                sim_times = simulated_values.get(f'{name}_times', None)
                
                ax.plot(target_times, target_values, 'b-', linewidth=2, label='Target', alpha=0.7)
                if sim_times is not None:
                    ax.plot(sim_times, sim_value, 'r--', linewidth=2, label='Simulated', alpha=0.7)
                else:
                    ax.plot(range(len(sim_value)), sim_value, 'r--', linewidth=2, label='Simulated', alpha=0.7)
                ax.set_xlabel('Time', fontsize=10)
                ax.set_ylabel(name, fontsize=10)
        else:
            # Plot scalar comparison
            target_value = target.get('target_value', 0.0)
            # Ensure sim_value is a scalar
            if isinstance(sim_value, np.ndarray):
                sim_value = float(sim_value.item() if sim_value.size == 1 else sim_value[0])
            else:
                sim_value = float(sim_value)
            target_value = float(target_value)
            ax.bar(['Target', 'Simulated'], [target_value, sim_value], 
                  color=['blue', 'red'], alpha=0.7)
            ax.set_ylabel('Value', fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.set_title(f'{name} ({target_type})', fontsize=11, fontweight='bold')
            continue  # Skip the duplicate title setting below
        
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Target: {name}', fontsize=11, fontweight='bold')
    
    # Hide unused subplots
    for idx in range(n_targets, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved target comparison plot to {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_simulation_results(
    results_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    output_dir: str,
    title_prefix: str = ""
):
    """
    Create comprehensive plots of simulation results.
    
    Args:
        results_df: DataFrame with time series results (must have 'time' column)
        summary_df: DataFrame with summary statistics (output_name, min, max, mean, std)
        output_dir: Directory to save plots
        title_prefix: Prefix for plot titles (e.g., "Baseline", "Optimized")
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Group outputs by category for better visualization
    pressure_outputs = [col for col in results_df.columns if 'pressure' in col.lower() and col != 'time']
    flow_outputs = [col for col in results_df.columns if 'flow' in col.lower() and col != 'time']
    volume_outputs = [col for col in results_df.columns if 
                      (col.lower().startswith('vc:') or col.lower().startswith('volume:') or 'volume' in col.lower()) 
                      and 'pressure' not in col.lower() and 'flow' not in col.lower() and col != 'time']
    
    # Plot 1: All Pressures
    if pressure_outputs:
        fig, ax = plt.subplots(figsize=(14, 8))
        for output in pressure_outputs:
            ax.plot(results_df['time'], results_df[output], label=output, linewidth=1.5)
        ax.set_xlabel('Time', fontsize=12)
        ax.set_ylabel('Pressure', fontsize=12)
        title = f'{title_prefix} Pressures ({len(pressure_outputs)} outputs)' if title_prefix else f'All Pressures ({len(pressure_outputs)} outputs)'
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'pressures.png'), dpi=150, bbox_inches='tight')
        plt.close()
    
    # Plot 2: All Flows
    if flow_outputs:
        fig, ax = plt.subplots(figsize=(14, 8))
        for output in flow_outputs:
            ax.plot(results_df['time'], results_df[output], label=output, linewidth=1.5)
        ax.set_xlabel('Time', fontsize=12)
        ax.set_ylabel('Flow', fontsize=12)
        title = f'{title_prefix} Flows ({len(flow_outputs)} outputs)' if title_prefix else f'All Flows ({len(flow_outputs)} outputs)'
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'flows.png'), dpi=150, bbox_inches='tight')
        plt.close()
    
    # Plot 3: All Volumes
    if volume_outputs:
        fig, ax = plt.subplots(figsize=(14, 8))
        for output in volume_outputs:
            ax.plot(results_df['time'], results_df[output], label=output, linewidth=1.5)
        ax.set_xlabel('Time', fontsize=12)
        ax.set_ylabel('Volume', fontsize=12)
        title = f'{title_prefix} Volumes ({len(volume_outputs)} outputs)' if title_prefix else f'All Volumes ({len(volume_outputs)} outputs)'
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'volumes.png'), dpi=150, bbox_inches='tight')
        plt.close()
    
    print(f"Saved simulation plots to: {output_dir}/")


def create_optimization_report(
    history: List[Dict],
    param_names: List[str],
    best_value: float,
    best_params: np.ndarray,
    output_dir: str
):
    """
    Create comprehensive optimization report with all visualizations.
    
    Args:
        history: Optimization history
        param_names: Parameter names
        best_value: Best objective value
        best_params: Best parameter values
        output_dir: Output directory for plots
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot objective history
    plot_objective_history(
        history,
        output_file=os.path.join(output_dir, 'objective_history.png')
    )
    
    # Plot parameter evolution
    if param_names:
        plot_parameter_evolution(
            history,
            param_names,
            output_file=os.path.join(output_dir, 'parameter_evolution.png')
        )
    
    print(f"\nOptimization report saved to {output_dir}")
    print(f"Best objective value: {best_value:.6e}")
    print(f"Best parameters:")
    for name, value in zip(param_names, best_params):
        print(f"  {name}: {value:.6e}")
