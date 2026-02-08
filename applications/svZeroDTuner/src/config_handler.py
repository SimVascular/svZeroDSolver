"""
Configuration Handler for sv0D Tuning Framework

Handles parsing and validation of YAML configuration files.
"""

import os
from typing import Dict, List, Any, Optional
from pathlib import Path

try:
    import yaml
except ImportError:
    raise ImportError(
        "PyYAML is required but not installed. Please install it with: pip install pyyaml"
    )


class ConfigHandler:
    """
    Handles loading and validation of YAML configuration files.
    """
    
    def __init__(self, config_file: str):
        """
        Initialize config handler.
        
        Args:
            config_file: Path to YAML configuration file
        """
        self.config_file = config_file
        self.config = self._load_config()
        self._validate_config()
    
    def _load_config(self) -> Dict:
        """Load YAML configuration file."""
        with open(self.config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Resolve relative paths relative to config file directory
        config_dir = os.path.dirname(os.path.abspath(self.config_file))
        if 'model' in config and 'config_file' in config['model']:
            if not os.path.isabs(config['model']['config_file']):
                config['model']['config_file'] = os.path.join(
                    config_dir, config['model']['config_file']
                )
        
        # Resolve target file paths
        if 'targets' in config:
            for target in config['targets']:
                if 'target_file' in target:
                    if not os.path.isabs(target['target_file']):
                        target['target_file'] = os.path.join(
                            config_dir, target['target_file']
                        )
        
        return config
    
    def _validate_config(self):
        """Validate configuration structure."""
        required_sections = ['model', 'parameters', 'targets', 'optimization']
        for section in required_sections:
            if section not in self.config:
                raise ValueError(f"Missing required section '{section}' in configuration")
        
        # Validate model section
        if 'config_file' not in self.config['model']:
            raise ValueError("model.config_file is required")
        if not os.path.exists(self.config['model']['config_file']):
            raise ValueError(f"Model config file not found: {self.config['model']['config_file']}")
        
        # Validate parameters
        if not isinstance(self.config['parameters'], list):
            raise ValueError("parameters must be a list")
        for param in self.config['parameters']:
            if 'name' not in param:
                raise ValueError("Each parameter must have 'name'")
            if 'bounds' not in param:
                raise ValueError(f"Parameter '{param['name']}' must have 'bounds'")
            if not isinstance(param['bounds'], list) or len(param['bounds']) != 2:
                raise ValueError(f"Parameter '{param['name']}' bounds must be [min, max]")
        
        # Validate targets
        if not isinstance(self.config['targets'], list):
            raise ValueError("targets must be a list")
        for target in self.config['targets']:
            if 'name' not in target:
                raise ValueError("Each target must have 'name'")
            if 'type' not in target:
                raise ValueError(f"Target '{target['name']}' must have 'type'")
            target_type = target['type']
            if target_type == 'time_series':
                if 'target_file' not in target:
                    raise ValueError(f"Time series target '{target['name']}' must have 'target_file'")
            elif target_type in ['min', 'max', 'mean']:
                if 'target_value' not in target:
                    raise ValueError(f"Scalar target '{target['name']}' must have 'target_value'")
            else:
                raise ValueError(f"Unknown target type: {target_type}")
        # Validate optimization section
        if 'algorithm' not in self.config['optimization']:
            raise ValueError("optimization.algorithm is required")
    
    def get_model_config_file(self) -> str:
        """Get path to sv0D.json model configuration file."""
        return self.config['model']['config_file']
    
    def get_parameters(self) -> List[Dict]:
        """Get list of parameters to optimize."""
        return self.config['parameters']
    
    def get_targets(self) -> List[Dict]:
        """Get list of targets."""
        return self.config['targets']
    
    def get_objective_config(self) -> Dict:
        """Get objective function configuration."""
        return self.config.get('objective', {'normalize': False})
    
    def get_optimization_config(self) -> Dict:
        """Get optimization configuration."""
        opt_config = self.config['optimization']
        return {
            'algorithm': opt_config.get('algorithm', 'differential_evolution'),
            'max_iterations': opt_config.get('max_iterations', 100),
            'tolerance': opt_config.get('tolerance', 1e-6),
            'parallel': opt_config.get('parallel', False),
            'n_workers': opt_config.get('n_workers', -1),
            **{k: v for k, v in opt_config.items() 
               if k not in ['algorithm', 'max_iterations', 'tolerance', 'parallel', 'n_workers']}
        }
    
    def get_output_config(self) -> Dict:
        """Get output configuration."""
        return self.config.get('output', {
            'directory': 'optimization_results',
            'save_history': True,
            'save_plots': True,
            'save_final_config': True
        })
    
    def get_config(self) -> Dict:
        """Get full configuration dictionary."""
        return self.config
