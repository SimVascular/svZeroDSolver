"""
Visualization functions for sv0D Tuning Framework

Plots optimization history including objective value convergence and parameter evolution.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
from scipy.interpolate import interp1d
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
    
    steps = [h.get('iteration', i) for i, h in enumerate(history)]
    objectives = [h['objective'] for h in history]
    
    plt.figure(figsize=(10, 6))
    plt.plot(steps, objectives, 'ko-', linewidth=2, markersize=6, label='Objective Value')
    plt.xlabel('Iteration', fontsize=12)
    plt.ylabel('Objective Value', fontsize=12)
    plt.title('Optimization Convergence', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    
    # Add best value annotation
    best_idx = np.argmin(objectives)
    best_step = steps[best_idx]
    best_obj = objectives[best_idx]
    plt.plot(best_step, best_obj, 'ro', markersize=10, label=f'Best: {best_obj:.6e}')
    plt.annotate(
        f'Best: {best_obj:.6e}',
        xy=(best_step, best_obj),
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
    
    steps = [h.get('iteration', i) for i, h in enumerate(history)]
    
    for idx, param_name in enumerate(param_names):
        ax = axes[idx]
        param_values = [h['parameters'][param_name] for h in history]
        
        ax.plot(steps, param_values, 'ko-', linewidth=2, markersize=6)
        ax.set_xlabel('Iteration', fontsize=10)
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


def _compute_percent_error(target_value: float, sim_value: float) -> Optional[float]:
    """Compute percent error: 100 * (sim - target) / target. Returns None if target ~ 0."""
    if abs(target_value) < 1e-14:
        return None
    return 100.0 * (sim_value - target_value) / target_value


def _get_range_for_target(target: Dict):
    """Get (lo, hi) for target. Returns None if not available."""
    if 'range_lo' in target and 'range_hi' in target:
        lo, hi = np.asarray(target['range_lo']), np.asarray(target['range_hi'])
        return (lo, hi)
    return None


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

    csv_rows = []
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
        
        if name not in simulated_values:
            ax.text(0.5, 0.5, f'No data for {name}', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Target: {name}', fontsize=11, fontweight='bold')
            csv_rows.append({'name': name, 'type': target_type, 'time': 'N/A', 'target_value': '', 'simulated_value': '', 'target_range': '', 'percent_error': 'N/A'})
            continue
        
        sim_value = simulated_values[name]
        
        if target_type == 'time_series':
            # Plot time series: use only range info; target_values = pointwise (lo+hi)/2
            rng = _get_range_for_target(target)
            if rng is None or 'target_times' not in target:
                ax.text(0.5, 0.5, f'No range data for {name}', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'Target: {name}', fontsize=11, fontweight='bold')
                csv_rows.append({'name': name, 'type': 'time_series', 'time': 'N/A', 'target_value': '', 'simulated_value': '', 'target_range': '', 'percent_error': 'N/A'})
                continue
            lo, hi = rng
            target_times = np.array(target['target_times'])
            target_values = (lo + hi) / 2.0  # pointwise midpoint
            sim_times = simulated_values.get(f'{name}_times', None)
            sim_val_arr = np.asarray(sim_value)
            ax.fill_between(target_times, lo, hi, alpha=0.2, color='blue', label='Target range')
            ax.plot(target_times, target_values, 'b-o', linewidth=2, markersize=4, label='Target', alpha=0.7)
            if sim_times is not None:
                ax.plot(sim_times, sim_value, 'r--', linewidth=2, label='Simulated', alpha=0.7)
            else:
                ax.plot(range(len(sim_val_arr)), sim_val_arr, 'r--', linewidth=2, label='Simulated', alpha=0.7)
            ax.set_xlabel('Time', fontsize=10)
            ax.set_ylabel(name, fontsize=10)
            # Compute MAPE and CSV rows (time series: one row per time point)
            if len(target_times) > 0 and len(sim_val_arr) > 1:
                try:
                    st = np.asarray(sim_times) if sim_times is not None else np.linspace(0, 1, len(sim_val_arr))
                    interp_func = interp1d(
                        st, sim_val_arr,
                        kind='linear', bounds_error=False, fill_value='extrapolate'
                    )
                    sim_interp = interp_func(target_times)
                    valid = np.isfinite(target_values) & (np.abs(target_values) > 1e-14)
                    if np.any(valid):
                        pct_errors = np.array([_compute_percent_error(float(target_values[i]), float(sim_interp[i])) for i in range(len(target_times))])
                        mape_vals = [abs(p) for p in pct_errors if p is not None]
                        mape = np.mean(mape_vals) if mape_vals else None
                        err_str = f'MAPE: {mape:.2f}%' if mape is not None else None
                        for i in range(len(target_times)):
                            range_str = f"[{float(lo[i]):.6e}, {float(hi[i]):.6e}]"
                            pct = pct_errors[i] if pct_errors[i] is not None else 'N/A'
                            csv_rows.append({
                                'name': name, 'type': 'time_series',
                                'time': target_times[i], 'target_value': target_values[i],
                                'simulated_value': sim_interp[i], 'target_range': range_str, 'percent_error': pct
                            })
                    else:
                        err_str = None
                except Exception:
                    err_str = None
            else:
                err_str = None
            if err_str is None and len(target_times) > 0 and len(sim_val_arr) > 0:
                st = np.asarray(sim_times) if sim_times is not None else np.linspace(0, 1, len(sim_val_arr))
                try:
                    interp_func = interp1d(st, sim_val_arr, kind='linear', bounds_error=False, fill_value='extrapolate')
                    sim_interp = interp_func(target_times)
                except Exception:
                    sim_interp = np.full_like(target_values, np.nan)
                pct_errors = [_compute_percent_error(float(target_values[i]), float(sim_interp[i])) for i in range(len(target_times))]
                for i in range(len(target_times)):
                    range_str = f"[{float(lo[i]):.6e}, {float(hi[i]):.6e}]"
                    pct = pct_errors[i] if pct_errors[i] is not None else 'N/A'
                    csv_rows.append({
                        'name': name, 'type': 'time_series',
                        'time': target_times[i], 'target_value': target_values[i],
                        'simulated_value': sim_interp[i], 'target_range': range_str, 'percent_error': pct
                    })
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
            ax.set_title(f'Target: {name}', fontsize=11, fontweight='bold')
            if err_str:
                ax.text(0.98, 0.05, f'Mean Absolute Percent Error: {mape:.2f}%', transform=ax.transAxes, fontsize=10,
                        ha='right', va='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            continue
        else:
            # Plot scalar: use only range info; target_value = (lo+hi)/2, error bar shows range
            rng = _get_range_for_target(target)
            if rng is None:
                ax.text(0.5, 0.5, f'No range data for {name}', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'Target: {name}', fontsize=11, fontweight='bold')
                csv_rows.append({'name': name, 'type': target_type, 'time': 'N/A', 'target_value': '', 'simulated_value': '', 'target_range': '', 'percent_error': 'N/A'})
                continue
            lo, hi = float(rng[0].flat[0]), float(rng[1].flat[0])
            target_value = (lo + hi) / 2.0
            range_str = f"[{lo:.6e}, {hi:.6e}]"
            if isinstance(sim_value, np.ndarray):
                sim_scalar = float(sim_value.item() if sim_value.size == 1 else sim_value[0])
            else:
                sim_scalar = float(sim_value)
            ax.bar(['Target', 'Simulated'], [target_value, sim_scalar], color=['blue', 'red'], alpha=0.7,
                   yerr=[[target_value - lo, 0], [hi - target_value, 0]], capsize=5)
            ax.set_ylabel('Value', fontsize=10)
            ax.grid(True, alpha=0.3)
            pct_err = _compute_percent_error(target_value, sim_scalar)
            err_str = f'Error: {pct_err:.2f}%' if pct_err is not None else None
            csv_rows.append({'name': name, 'type': target_type, 'time': 'N/A', 'target_value': target_value, 'simulated_value': sim_scalar, 'target_range': range_str, 'percent_error': pct_err if pct_err is not None else 'N/A'})
            ax.set_title(name, fontsize=11, fontweight='bold')
            if err_str:
                trans = blended_transform_factory(ax.transData, ax.transAxes)
                ax.text(1, 0.05, err_str, transform=trans, fontsize=10,
                       ha='center', va='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            continue
    
    # Hide unused subplots
    for idx in range(n_targets, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved target comparison plot to {output_file}")
        if csv_rows:
            csv_path = os.path.splitext(output_file)[0] + '.csv'
            col_order = ['name', 'type', 'time', 'target_value', 'simulated_value', 'target_range', 'percent_error']
            pd.DataFrame(csv_rows).to_csv(csv_path, index=False, columns=col_order)
            print(f"Saved target comparison data to {csv_path}")
    
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
    pressure_outputs = [col for col in results_df.columns
                        if (col.lower().startswith('pressure') or col.lower().startswith('pressure_c'))
                        and col != 'time']
    flow_outputs = [col for col in results_df.columns
                    if col.lower().startswith('flow:')
                    and col != 'time']
    volume_outputs = [col for col in results_df.columns if 
                      (col.lower().startswith('vc:') or col.lower().startswith('volume:') or 'volume' in col.lower()) 
                      and 'pressure' not in col.lower() and 'flow' not in col.lower() and col != 'time']
    
    # Use default axes color cycle; cycle linestyle when colors repeat
    default_colors = plt.rcParams['axes.prop_cycle'].by_key().get('color', list(plt.cm.tab10.colors))
    linestyles = ['-', '--', ':']

    def _plot_outputs(ax, outputs):
        for i, output in enumerate(outputs):
            color = default_colors[i % len(default_colors)]
            ls = linestyles[(i // len(default_colors)) % len(linestyles)]
            ax.plot(results_df['time'], results_df[output], label=output, linewidth=1.5, color=color, linestyle=ls)

    # Plot 1: All Pressures
    if pressure_outputs:
        fig, ax = plt.subplots(figsize=(14, 8))
        _plot_outputs(ax, pressure_outputs)
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
        _plot_outputs(ax, flow_outputs)
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
        _plot_outputs(ax, volume_outputs)
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
