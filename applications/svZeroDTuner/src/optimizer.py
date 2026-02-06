"""
Optimizer wrappers for sv0D Tuning Framework

Supports scipy optimizers and parallelization.
"""

import numpy as np
from typing import Callable, List, Tuple, Dict, Optional, Any
from scipy.optimize import minimize, differential_evolution, OptimizeResult
import multiprocessing as mp
from functools import partial


class OptimizerWrapper:
    """
    Wrapper around optimization algorithms with support for parallelization.
    """
    
    # Supported algorithms and their properties
    SUPPORTED_ALGORITHMS = {
        'differential_evolution': {
            'supports_bounds': True,
            'requires_x0': False,
            'supports_parallel': True,
            'valid_options': ['strategy', 'popsize', 'mutation', 'recombination', 'seed', 'polish', 'init']
        },
        'Nelder-Mead': {
            'supports_bounds': False,  # Enforced via penalty
            'requires_x0': True,
            'supports_parallel': False,
            'valid_options': ['initial_simplex', 'adaptive']
        },
        'L-BFGS-B': {
            'supports_bounds': True,
            'requires_x0': True,
            'supports_parallel': False,
            'valid_options': ['maxcor', 'maxls']
        },
        'BFGS': {
            'supports_bounds': False,
            'requires_x0': True,
            'supports_parallel': False,
            'valid_options': ['norm']
        },
        'Powell': {
            'supports_bounds': False,
            'requires_x0': True,
            'supports_parallel': False,
            'valid_options': ['direc']
        },
        'CG': {
            'supports_bounds': False,
            'requires_x0': True,
            'supports_parallel': False,
            'valid_options': []
        }
    }
    
    def __init__(
        self,
        algorithm: str = "differential_evolution",
        max_iterations: int = 100,
        tolerance: float = 1e-6,
        parallel: bool = False,
        n_workers: int = -1,
        **algorithm_kwargs
    ):
        """
        Initialize optimizer wrapper.
        
        Args:
            algorithm: Optimization algorithm name
            max_iterations: Maximum number of iterations
            tolerance: Convergence tolerance
            parallel: Whether to enable parallel evaluation
            n_workers: Number of parallel workers (-1 for all cores)
            **algorithm_kwargs: Additional algorithm-specific arguments
        """
        # Ensure proper types (in case values come from YAML as strings)
        self.algorithm = str(algorithm)
        self.max_iterations = int(max_iterations)
        self.tolerance = float(tolerance)
        self.parallel = bool(parallel)
        self.n_workers = int(n_workers) if n_workers > 0 else mp.cpu_count()
        self.algorithm_kwargs = algorithm_kwargs
        
        # Validate configuration
        self._validate_config()
        
        # History tracking
        self.history = []
        self.best_value = None
        self.best_params = None
        self.bounds = None  # Store bounds for enforcement
    
    def _validate_config(self):
        """
        Validate optimizer configuration.
        
        Raises:
            ValueError: If configuration is invalid
        """
        # Check if algorithm is known
        if self.algorithm not in self.SUPPORTED_ALGORITHMS:
            supported = ', '.join(self.SUPPORTED_ALGORITHMS.keys())
            print(f"Warning: Algorithm '{self.algorithm}' is not in the list of validated algorithms.")
            print(f"Supported algorithms: {supported}")
            print("Attempting to use it anyway, but this may fail.")
            return
        
        algo_config = self.SUPPORTED_ALGORITHMS[self.algorithm]
        
        # Check parallel support
        if self.parallel and not algo_config['supports_parallel']:
            raise ValueError(
                f"Algorithm '{self.algorithm}' does not support parallel evaluation. "
                f"Please set 'parallel: false' in tuning.yaml or choose a different algorithm "
                f"(e.g., 'differential_evolution')."
            )
        
        # Warn about bounds support
        if not algo_config['supports_bounds']:
            if self.algorithm == 'Nelder-Mead':
                print(f"Note: {self.algorithm} does not natively support bounds.")
                print("Bounds will be enforced via penalty function.")
            else:
                print(f"Warning: {self.algorithm} does not support bounds.")
                print("Parameters may go outside specified bounds during optimization.")
        
        # Check for invalid algorithm-specific options
        if self.algorithm_kwargs:
            valid_opts = algo_config['valid_options']
            invalid_opts = [k for k in self.algorithm_kwargs.keys() if k not in valid_opts]
            if invalid_opts and valid_opts:  # Only warn if we have a list of valid options
                print(f"Warning: Unknown options for {self.algorithm}: {invalid_opts}")
                if valid_opts:
                    print(f"Valid options are: {valid_opts}")
        
        # Validate parameter values
        if self.max_iterations <= 0:
            raise ValueError(f"max_iterations must be positive, got {self.max_iterations}")
        
        if self.tolerance <= 0:
            raise ValueError(f"tolerance must be positive, got {self.tolerance}")
        
        if self.parallel and self.n_workers <= 0:
            raise ValueError(f"n_workers must be positive when parallel=True, got {self.n_workers}")
    
    def _validate_optimization_inputs(
        self,
        bounds: List[Tuple[float, float]],
        x0: Optional[np.ndarray]
    ):
        """
        Validate optimization problem inputs.
        
        Args:
            bounds: Parameter bounds
            x0: Initial guess
            
        Raises:
            ValueError: If inputs are invalid
        """
        # Check bounds if provided
        if bounds:
            for i, (lower, upper) in enumerate(bounds):
                if lower >= upper:
                    raise ValueError(
                        f"Invalid bounds for parameter {i}: lower bound ({lower}) "
                        f"must be less than upper bound ({upper})"
                    )
        
        # Check if x0 is required and provided
        if self.algorithm in self.SUPPORTED_ALGORITHMS:
            algo_config = self.SUPPORTED_ALGORITHMS[self.algorithm]
            if algo_config['requires_x0'] and x0 is None:
                print(f"Note: {self.algorithm} typically requires an initial guess (x0).")
                print(f"Using center of bounds as initial guess.")
        
        # If x0 is provided, validate it
        if x0 is not None:
            if len(x0) != len(bounds):
                raise ValueError(
                    f"Initial guess x0 has {len(x0)} elements, "
                    f"but {len(bounds)} parameters specified"
                )
            
            # Check if x0 is within bounds
            for i, (x_val, (lower, upper)) in enumerate(zip(x0, bounds)):
                if x_val < lower or x_val > upper:
                    print(
                        f"Warning: Initial guess for parameter {i} ({x_val}) "
                        f"is outside bounds [{lower}, {upper}]. "
                        f"This may cause issues with some optimizers."
                    )
    
    def _objective_wrapper(
        self,
        params: np.ndarray,
        objective_func: Callable,
        param_names: List[str],
        iteration_callback: Optional[Callable] = None
    ) -> float:
        """
        Wrapper for objective function that tracks history and enforces bounds.
        
        Args:
            params: Parameter values
            objective_func: Objective function to call
            param_names: Names of parameters
            iteration_callback: Optional callback function for each iteration
            
        Returns:
            Objective function value (always a float)
        """
        # Check bounds for algorithms that don't natively support them
        if self.bounds is not None:
            out_of_bounds = False
            for i, (param, (lower, upper)) in enumerate(zip(params, self.bounds)):
                if param < lower or param > upper:
                    out_of_bounds = True
                    break
            
            if out_of_bounds:
                # Return large penalty value for out-of-bounds parameters
                return 1e10
        
        obj_value = objective_func(params)
        
        # Ensure obj_value is a float (not numpy scalar or array)
        obj_value = float(obj_value)
        
        # Track history
        history_entry = {
            'iteration': len(self.history),
            'objective': obj_value,
            'parameters': dict(zip(param_names, params.tolist()))  # Convert to list for JSON serialization
        }
        self.history.append(history_entry)
        
        # Update best
        if self.best_value is None or obj_value < self.best_value:
            self.best_value = obj_value
            self.best_params = params.copy()
        
        # Call iteration callback if provided
        if iteration_callback:
            iteration_callback(history_entry)
        
        return obj_value
    
    def optimize(
        self,
        objective_func: Callable,
        param_names: List[str],
        bounds: List[Tuple[float, float]],
        x0: Optional[np.ndarray] = None,
        iteration_callback: Optional[Callable] = None
    ) -> OptimizeResult:
        """
        Run optimization.
        
        Args:
            objective_func: Objective function that takes parameter array and returns scalar
            param_names: List of parameter names
            bounds: List of (min, max) tuples for each parameter
            x0: Initial guess (optional, required for some algorithms)
            iteration_callback: Optional callback function for each iteration
            
        Returns:
            Optimization result from scipy.optimize
        """
        # Reset history
        self.history = []
        self.best_value = None
        self.best_params = None
        
        # Store bounds for enforcement in objective wrapper
        self.bounds = bounds
        
        # Validate optimization inputs
        self._validate_optimization_inputs(bounds, x0)
        
        # Create wrapped objective
        wrapped_obj = partial(
            self._objective_wrapper,
            objective_func=objective_func,
            param_names=param_names,
            iteration_callback=iteration_callback
        )
        
        # Convert bounds to scipy format
        scipy_bounds = bounds
        
        if self.algorithm == "differential_evolution":
            result = differential_evolution(
                wrapped_obj,
                bounds=scipy_bounds,
                maxiter=self.max_iterations,
                tol=self.tolerance,
                workers=self.n_workers if self.parallel else 1,
                **self.algorithm_kwargs
            )
        
        elif self.algorithm == "Nelder-Mead":
            if x0 is None:
                # Use center of bounds as initial guess
                x0 = np.array([(b[0] + b[1]) / 2 for b in bounds])
            # Ensure tolerance is a float
            tol = float(self.tolerance)
            # Nelder-Mead method options (initial_simplex, adaptive) must go in options dict, not as top-level kwargs
            nm_options = {
                'maxiter': int(self.max_iterations),
                'xatol': tol,
                'fatol': tol
            }
            valid_nm_opts = self.SUPPORTED_ALGORITHMS['Nelder-Mead']['valid_options']
            for k, v in self.algorithm_kwargs.items():
                if k in valid_nm_opts:
                    nm_options[k] = v
            # Note: Nelder-Mead doesn't natively support bounds, so we enforce them via penalty in objective wrapper
            result = minimize(
                wrapped_obj,
                x0=x0,
                method='Nelder-Mead',
                options=nm_options
            )
        
        elif self.algorithm == "L-BFGS-B":
            if x0 is None:
                x0 = np.array([(b[0] + b[1]) / 2 for b in bounds])
            result = minimize(
                wrapped_obj,
                x0=x0,
                method='L-BFGS-B',
                bounds=scipy_bounds,
                options={
                    'maxiter': self.max_iterations,
                    'gtol': self.tolerance
                },
                **self.algorithm_kwargs
            )
        
        elif self.algorithm == "BFGS":
            if x0 is None:
                x0 = np.array([(b[0] + b[1]) / 2 for b in bounds])
            result = minimize(
                wrapped_obj,
                x0=x0,
                method='BFGS',
                options={
                    'maxiter': self.max_iterations,
                    'gtol': self.tolerance
                },
                **self.algorithm_kwargs
            )
        
        else:
            # Try to use as scipy method name directly
            if x0 is None:
                x0 = np.array([(b[0] + b[1]) / 2 for b in bounds])
            result = minimize(
                wrapped_obj,
                x0=x0,
                method=self.algorithm,
                bounds=scipy_bounds if scipy_bounds else None,
                options={
                    'maxiter': self.max_iterations
                },
                **self.algorithm_kwargs
            )
        
        return result
    
    def get_history(self) -> List[Dict]:
        """Get optimization history."""
        return self.history
    
    def get_best(self) -> Tuple[float, np.ndarray]:
        """Get best objective value and parameters."""
        return self.best_value, self.best_params
    
    @classmethod
    def print_supported_algorithms(cls):
        """Print information about supported algorithms."""
        print("\n" + "="*70)
        print("SUPPORTED OPTIMIZATION ALGORITHMS")
        print("="*70)
        for algo, config in cls.SUPPORTED_ALGORITHMS.items():
            print(f"\n{algo}:")
            print(f"  - Supports bounds: {config['supports_bounds']}")
            print(f"  - Requires initial guess: {config['requires_x0']}")
            print(f"  - Supports parallel: {config['supports_parallel']}")
            if config['valid_options']:
                print(f"  - Valid options: {', '.join(config['valid_options'])}")
        print("\n" + "="*70)
