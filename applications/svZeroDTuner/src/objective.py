"""
Objective Function implementations for sv0D Tuning Framework

Supports multiple objective function types: L2 norm, relative error, weighted combinations.
"""

import numpy as np
from typing import Dict, List, Callable, Optional, Union
from scipy.interpolate import interp1d


class ObjectiveFunction:
    """
    Base class for objective functions.
    """
    
    def __init__(
        self,
        targets: List[Dict],
        normalize: bool = False,
        custom_function: Optional[Callable] = None
    ):
        """
        Initialize objective function.
        
        Args:
            targets: List of target specifications, each with:
                - name: Output name
                - type: Extraction type (time_series, min, max, mean)
                - target_value or target_file: Target value(s)
                - weight: Weight for this target (default 1.0)
            normalize: Whether to normalize errors by target values
            custom_function: Custom objective function (if provided, overrides default)
        """
        self.targets = targets
        self.normalize = normalize
        self.custom_function = custom_function
        
        # Process targets
        self._process_targets()
    
    def _process_targets(self):
        """Process and validate target specifications."""
        for target in self.targets:
            if 'weight' not in target:
                target['weight'] = 1.0
            
            # Load target values
            if 'target_file' in target:
                # Load time series from CSV
                import pandas as pd
                df = pd.read_csv(target['target_file'])
                if 'time' not in df.columns or 'value' not in df.columns:
                    raise ValueError(f"target_file must have 'time' and 'value' columns: {target['target_file']}")
                target['target_times'] = df['time'].values
                target['target_values'] = df['value'].values
            elif 'target_value' in target:
                target['target_value'] = float(target['target_value'])
            else:
                raise ValueError(f"Target must have either 'target_file' or 'target_value': {target}")
    
    def compute_error(
        self,
        simulated_values: Dict[str, Union[np.ndarray, float]]
    ) -> float:
        """
        Compute objective function value.
        
        Args:
            simulated_values: Dictionary mapping output names to simulated values
            
        Returns:
            Objective function value (to minimize)
        """
        if self.custom_function:
            return self.custom_function(simulated_values, self.targets)
        
        total_error = 0.0
        
        for target in self.targets:
            name = target['name']
            extraction_type = target.get('type', 'time_series')
            weight = float(target.get('weight', 1.0))
            
            # Create unique key for this target (name + type to handle duplicates)
            target_key = f"{name}_{extraction_type}" if extraction_type != 'time_series' else name
            
            # For time_series, use the base name; for others, try the composite key or base name
            if target_key in simulated_values:
                sim_value = simulated_values[target_key]
            elif name in simulated_values:
                sim_value = simulated_values[name]
            else:
                raise ValueError(f"Simulated value for '{name}' (type: {extraction_type}) not found")
            
            # Ensure sim_value is numpy array
            if not isinstance(sim_value, np.ndarray):
                sim_value = np.array([sim_value]) if np.isscalar(sim_value) else np.array(sim_value)
            
            if extraction_type == 'time_series':
                # Compare time series
                if 'target_times' not in target:
                    raise ValueError(f"Time series target '{name}' must have 'target_file'")
                
                target_times = target['target_times']
                target_values = target['target_values']
                
                # Interpolate simulated values to target time points
                # Try to get times from simulated_values, or use indices
                sim_times = simulated_values.get(f'{name}_times', None)
                if sim_times is None:
                    # If no times provided, assume uniform spacing
                    sim_times = np.linspace(0, 1, len(sim_value))
                
                # Ensure arrays are numpy arrays and have same length
                sim_times = np.array(sim_times)
                sim_value = np.array(sim_value)
                target_times = np.array(target_times)
                target_values = np.array(target_values)
                
                # Check array lengths match
                if len(sim_times) != len(sim_value):
                    raise ValueError(
                        f"Time and value arrays have different lengths for {name}: "
                        f"{len(sim_times)} vs {len(sim_value)}"
                    )
                
                # Interpolate
                if len(sim_times) > 1 and len(target_times) > 0:
                    # Remove any NaN or inf values
                    valid_mask = np.isfinite(sim_value) & np.isfinite(sim_times)
                    if np.sum(valid_mask) < 2:
                        # Not enough valid points, return large error
                        error = 1e10
                    else:
                        interp_func = interp1d(
                            sim_times[valid_mask], 
                            sim_value[valid_mask], 
                            kind='linear', 
                            bounds_error=False, 
                            fill_value='extrapolate'
                        )
                        sim_interp = interp_func(target_times)
                        # Compute error
                        if self.normalize:
                            # Normalized L2 error
                            norm = np.linalg.norm(target_values)
                            if norm > 0:
                                error = np.linalg.norm(sim_interp - target_values) / norm
                            else:
                                error = np.linalg.norm(sim_interp - target_values)
                        else:
                            # L2 error
                            error = np.linalg.norm(sim_interp - target_values) ** 2
                else:
                    # Not enough points for interpolation
                    error = 1e10
                
                # Error already computed above
                pass
            
            else:
                # Scalar comparison (min, max, mean)
                if 'target_value' not in target:
                    raise ValueError(f"Scalar target '{name}' must have 'target_value'")
                
                target_value = target['target_value']
                
                # Compute error
                if self.normalize:
                    # Relative error
                    if abs(target_value) > 1e-10:
                        error = abs(sim_value - target_value) / abs(target_value)
                    else:
                        error = abs(sim_value - target_value)
                else:
                    # Absolute error squared
                    error = (sim_value - target_value) ** 2
            
            total_error += weight * error
        
        # Ensure we return a float, not numpy scalar or array
        if isinstance(total_error, np.ndarray):
            return float(total_error.item())  # .item() will error if size > 1
        return float(total_error)


class L2Objective(ObjectiveFunction):
    """L2 norm objective function (sum of squared differences)."""
    
    def __init__(self, targets: List[Dict], **kwargs):
        super().__init__(targets, normalize=False, **kwargs)


class RelativeErrorObjective(ObjectiveFunction):
    """Relative error objective function (normalized by target values)."""
    
    def __init__(self, targets: List[Dict], **kwargs):
        super().__init__(targets, normalize=True, **kwargs)


class WeightedL2Objective(ObjectiveFunction):
    """Weighted L2 objective function (uses weights from target specifications)."""
    
    def __init__(self, targets: List[Dict], **kwargs):
        super().__init__(targets, **kwargs)


def create_objective(
    objective_type: str,
    targets: List[Dict],
    **kwargs
) -> ObjectiveFunction:
    """
    Factory function to create objective function.
    
    Args:
        objective_type: Type of objective ('l2', 'relative_error', 'weighted_l2', 'custom')
        targets: List of target specifications
        **kwargs: Additional arguments for objective function
        
    Returns:
        ObjectiveFunction instance
    """
    if objective_type == 'l2':
        return L2Objective(targets, **kwargs)
    elif objective_type == 'relative_error':
        return RelativeErrorObjective(targets, **kwargs)
    elif objective_type == 'weighted_l2':
        return WeightedL2Objective(targets, **kwargs)
    elif objective_type == 'custom':
        return ObjectiveFunction(targets, custom_function=kwargs.get('custom_function'), **kwargs)
    else:
        raise ValueError(f"Unknown objective type: {objective_type}")
