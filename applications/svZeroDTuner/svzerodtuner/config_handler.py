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
        required_sections = ['model', 'parameters', 'targets', 'objective', 'optimization']
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
            # Optional scaling: "identity", "log", or "max"
            scaling = param.get('scaling')
            if scaling is not None:
                if scaling == 'identity':
                    pass  # no extra validation
                elif scaling == 'log':
                    lo, hi = float(param['bounds'][0]), float(param['bounds'][1])
                    if lo <= 0 or hi <= 0:
                        raise ValueError(
                            f"Parameter '{param['name']}' has scaling 'log'; bounds must be positive, got [{lo}, {hi}]"
                        )
                elif scaling == 'max':
                    lo, hi = float(param['bounds'][0]), float(param['bounds'][1])
                    if max(lo, hi) == 0:
                        raise ValueError(
                            f"Parameter '{param['name']}' has scaling 'max'; bounds must be non-zero"
                        )
                else:
                    raise ValueError(
                        f"Parameter '{param['name']}' scaling must be 'identity', 'log', or 'max' (got {scaling!r})"
                    )
        
        # Validate targets
        if not isinstance(self.config['targets'], list):
            raise ValueError("targets must be a list")
        for target in self.config['targets']:
            if 'name' not in target:
                raise ValueError("Each target must have 'name'")
            if 'type' not in target:
                raise ValueError(f"Target '{target['name']}' must have 'type'")
            has_relative = 'relative_bounds' in target
            has_uncertainty = 'uncertainty' in target
            if has_relative and has_uncertainty:
                raise ValueError(
                    f"Target '{target['name']}' cannot define both 'relative_bounds' and legacy 'uncertainty'; use only one"
                )
            if has_relative or has_uncertainty:
                unc = target['relative_bounds'] if has_relative else target['uncertainty']
                if isinstance(unc, str) and unc.strip().endswith('%'):
                    raw = unc.strip()[:-1]
                    try:
                        pct = float(raw)
                    except ValueError:
                        raise ValueError(
                            f"Target '{target['name']}' relative_bounds '{unc}' must be a valid percent (e.g. '5%')"
                        )
                    if pct < 0:
                        raise ValueError(
                            f"Target '{target['name']}' relative_bounds percent must be non-negative"
                        )
                elif isinstance(unc, (int, float)):
                    if unc < 0:
                        raise ValueError(
                            f"Target '{target['name']}' relative_bounds percent must be non-negative"
                        )
                elif isinstance(unc, (list, tuple)):
                    key_name = 'relative_bounds' if has_relative else 'uncertainty'
                    if len(unc) != 2:
                        raise ValueError(
                            f"Target '{target['name']}' {key_name} [min, max] must have 2 elements"
                        )
                    if unc[0] >= unc[1]:
                        raise ValueError(
                            f"Target '{target['name']}' {key_name} [min, max] must have min < max"
                        )
                else:
                    raise ValueError(
                        f"Target '{target['name']}' relative_bounds must be percent (e.g. '5%') or [min, max]"
                    )
            target_type = target['type']
            if target_type == 'time_series':
                if 'expression' not in target:
                    raise ValueError(
                        f"Time series target '{target['name']}' must have 'expression'"
                    )
                if 'target_file' not in target:
                    raise ValueError(
                        f"Time series target '{target['name']}' must have 'target_file'"
                    )
                if 'target_range' in target:
                    r = target['target_range']
                    if not isinstance(r, (list, tuple)) or len(r) != 2:
                        raise ValueError(
                            f"Target '{target['name']}' target_range must be [min, max]"
                        )
                    if r[0] >= r[1]:
                        raise ValueError(
                            f"Target '{target['name']}' target_range must have min < max"
                        )
            elif target_type == 'scalar':
                if 'expression' not in target:
                    raise ValueError(
                        f"Scalar target '{target['name']}' must have 'expression'"
                    )
                if 'target_value' not in target and 'target_range' not in target:
                    raise ValueError(
                        f"Scalar target '{target['name']}' must have 'target_value' or 'target_range'"
                    )
                if 'target_range' in target:
                    r = target['target_range']
                    if not isinstance(r, (list, tuple)) or len(r) != 2:
                        raise ValueError(
                            f"Target '{target['name']}' target_range must be [min, max]"
                        )
                    if r[0] >= r[1]:
                        raise ValueError(
                            f"Target '{target['name']}' target_range must have min < max"
                        )
            else:
                raise ValueError(f"Unknown target type: {target_type}")
        # Validate optimization section
        if 'algorithm' not in self.config['optimization']:
            raise ValueError("optimization.algorithm is required")
        
        # Validate objective section (norm is required)
        if not isinstance(self.config.get('objective'), dict):
            raise ValueError(
                "objective section is required and must be a mapping. "
                "Example:\n  objective:\n    norm: L1\n"
                "Options for norm: L1 = sum of absolute relative errors; L2 = Euclidean norm of the error vector."
            )
        if 'norm' not in self.config['objective']:
            raise ValueError(
                "objective.norm is required"
            )
    
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
        """Get objective function configuration (e.g. norm: 'L1' or 'L2')."""
        return self.config['objective']
    
    def get_optimization_config(self) -> Dict:
        """Get optimization config. Passed directly to optimizer"""
        return dict(self.config['optimization'])
    
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
