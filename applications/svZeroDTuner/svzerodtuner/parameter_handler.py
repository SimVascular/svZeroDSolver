"""
Parameter Handler for sv0D Tuning Framework

Handles loading, parsing, and updating parameters in sv0D.json files using name-based lookup.
"""

import json
import copy
import numpy as np
from typing import Dict, List, Tuple, Any, Optional


class ParameterHandler:
    """
    Handles parameter loading, lookup, and updating in sv0D JSON configuration files.
    
    Supports name-based parameter lookup (e.g., "LV.Emax", "AR_SYS.C")
    for chambers, vessels, valves, and other structures.
    """
    
    def __init__(self, config_file: str):
        """
        Initialize parameter handler with a sv0D.json file.
        
        Args:
            config_file: Path to sv0D.json configuration file
        """
        self.config_file = config_file
        self.config = self._load_config()
        self.original_config = copy.deepcopy(self.config)
    
    def _load_config(self) -> Dict:
        """Load JSON configuration file."""
        with open(self.config_file, 'r') as f:
            return json.load(f)
    
    def find_parameter(self, param_name: str) -> Tuple[Any, List[Any]]:
        """
        Find a parameter by name using name-based lookup.
        
        Supports formats like:
        - "LV.Emax" -> chambers[name="LV"].values.Emax
        - "AR_SYS.C" -> vessels[name="AR_SYS"].zero_d_element_values.C
        - "MV.Rmin" -> valves[name="MV"].params.Rmin
        
        Args:
            param_name: Parameter name in format "BlockName.ParameterName"
            
        Returns:
            Tuple of (value, path) where path is list of keys/indexes to reach the parameter
            
        Raises:
            ValueError: If parameter is not found
        """
        parts = param_name.split('.')
        if len(parts) != 2:
            raise ValueError(f"Parameter name must be in format 'BlockName.ParameterName', got '{param_name}'")
        
        block_name, param_key = parts
        
        # Try chambers first
        if 'chambers' in self.config:
            for i, chamber in enumerate(self.config['chambers']):
                if chamber.get('name') == block_name:
                    if 'values' in chamber and param_key in chamber['values']:
                        return chamber['values'][param_key], ['chambers', i, 'values', param_key]
        
        # Try vessels
        if 'vessels' in self.config:
            for i, vessel in enumerate(self.config['vessels']):
                if vessel.get('vessel_name') == block_name:
                    # Check zero_d_element_values
                    if 'zero_d_element_values' in vessel:
                        if param_key in vessel['zero_d_element_values']:
                            return vessel['zero_d_element_values'][param_key], ['vessels', i, 'zero_d_element_values', param_key]
                    # Check other top-level vessel properties
                    if param_key in vessel:
                        return vessel[param_key], ['vessels', i, param_key]
        
        # Try valves
        if 'valves' in self.config:
            for i, valve in enumerate(self.config['valves']):
                if valve.get('name') == block_name:
                    if 'params' in valve and param_key in valve['params']:
                        return valve['params'][param_key], ['valves', i, 'params', param_key]
        
        # Try boundary conditions
        if 'boundary_conditions' in self.config:
            for i, bc in enumerate(self.config['boundary_conditions']):
                if bc.get('bc_name') == block_name:
                    if 'bc_values' in bc and param_key in bc['bc_values']:
                        return bc['bc_values'][param_key], ['boundary_conditions', i, 'bc_values', param_key]
        
        # Try simulation parameters
        if 'simulation_parameters' in self.config:
            if block_name == 'simulation' and param_key in self.config['simulation_parameters']:
                return self.config['simulation_parameters'][param_key], ['simulation_parameters', param_key]
        
        # Try initial conditions
        if 'initial_condition' in self.config:
            if block_name == 'initial_condition' and param_key in self.config['initial_condition']:
                return self.config['initial_condition'][param_key], ['initial_condition', param_key]
        
        raise ValueError(f"Parameter '{param_name}' not found in configuration")
    
    def set_parameter(self, param_name: str, value: Any) -> None:
        """
        Set a parameter value by name.
        
        Args:
            param_name: Parameter name in format "BlockName.ParameterName"
            value: New parameter value (will be converted to Python scalar if numpy array)
        """
        # Convert numpy arrays/types to Python scalars
        if isinstance(value, np.ndarray):
            value = float(value.item())
        elif isinstance(value, (np.integer, np.floating)):
            value = float(value)

        _, path = self.find_parameter(param_name)
        target = self.config
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
    
    def get_parameter(self, param_name: str) -> Any:
        """
        Get a parameter value by name.
        
        Args:
            param_name: Parameter name in format "BlockName.ParameterName"
            
        Returns:
            Parameter value
        """
        value, _ = self.find_parameter(param_name)
        return value
    
    def _convert_numpy_types(self, obj: Any) -> Any:
        """
        Recursively convert numpy types to Python native types.
        
        Args:
            obj: Object to convert
            
        Returns:
            Object with numpy types converted to Python types
        """
        if isinstance(obj, np.ndarray):
            return float(obj.item()) if obj.size == 1 else obj.tolist()
        elif isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        elif isinstance(obj, dict):
            return {key: self._convert_numpy_types(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [self._convert_numpy_types(item) for item in obj]
        else:
            return obj
    
    def get_config(self) -> Dict:
        """
        Get the current configuration dictionary with all numpy types converted to Python types.
        
        Returns:
            Configuration dictionary
        """
        config_copy = copy.deepcopy(self.config)
        return self._convert_numpy_types(config_copy)
    
    def reset_to_original(self) -> None:
        """Reset configuration to original loaded values."""
        self.config = copy.deepcopy(self.original_config)
