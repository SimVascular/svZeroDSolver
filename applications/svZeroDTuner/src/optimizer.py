"""
Optimizer wrappers for sv0D Tuning Framework

Supports scipy optimizers and parallelization.
"""

import numpy as np
from typing import Callable, List, Tuple, Dict, Optional, Any
from scipy.optimize import minimize, differential_evolution, OptimizeResult
from functools import partial


def _coerce_numeric_options(options: Dict) -> Dict:
    """Convert string values that look like numbers to int/float (YAML can load 1e-6, 1e-2 as str)."""
    def coerce(v):
        if isinstance(v, str):
            try:
                f = float(v)
                return int(f) if f.is_integer() else f
            except ValueError:
                return v
        if isinstance(v, (list, tuple)):
            return type(v)(coerce(x) for x in v)
        return v

    return {k: coerce(v) for k, v in options.items()}


class OptimizerWrapper:
    """
    Wrapper around optimization algorithms with support for parallelization.
    """
    
    SUPPORTED_ALGORITHMS = {'differential_evolution', 'Nelder-Mead'}
    
    def __init__(self, algorithm: str = "differential_evolution", **options):
        """
        Initialize optimizer wrapper.
        
        Args:
            algorithm: Optimization algorithm name ('differential_evolution' or 'Nelder-Mead')
            **options: All other options passed directly to the scipy optimization function.
                       Use scipy's exact parameter names. Invalid options will raise from scipy.
        """
        self.algorithm = str(algorithm)
        self.options = options
        self._validate_config()
        self.history = []
        self.best_value = None
        self.best_params = None
        self.bounds = None
    
    def _validate_config(self):
        if self.algorithm not in self.SUPPORTED_ALGORITHMS:
            raise ValueError(
                f"Algorithm '{self.algorithm}' is not supported. "
                f"Supported: {', '.join(sorted(self.SUPPORTED_ALGORITHMS))}"
            )
        if self.algorithm == 'Nelder-Mead':
            print("Note: Nelder-Mead does not natively support bounds; enforcing via penalty.")
    
    def _validate_optimization_inputs(
        self,
        bounds: List[Tuple[float, float]],
        x0: Optional[np.ndarray],
        param_names: Optional[List[str]] = None
    ):
        """
        Validate optimization problem inputs.
        
        Args:
            bounds: Parameter bounds
            x0: Initial guess
            param_names: Parameter names (optional, for clearer error messages)
            
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
        
        if self.algorithm == 'Nelder-Mead' and x0 is None:
            print("Note: Using center of bounds as initial guess for Nelder-Mead.")
        
        # If x0 is provided, validate it
        if x0 is not None:
            if len(x0) != len(bounds):
                raise ValueError(
                    f"Initial guess x0 has {len(x0)} elements, "
                    f"but {len(bounds)} parameters specified"
                )
            
            # Check if x0 is within bounds - raise error if any initial value is outside
            for i, (x_val, (lower, upper)) in enumerate(zip(x0, bounds)):
                if x_val < lower or x_val > upper:
                    name = param_names[i] if param_names and i < len(param_names) else f"parameter {i}"
                    bound_type = "min" if x_val < lower else "max"
                    bound_val = lower if x_val < lower else upper
                    raise ValueError(
                        f"Initial value for {name} ({x_val}) is outside its bounds: "
                        f"violates {bound_type} bound ({bound_val}). "
                        f"Bounds are [{lower}, {upper}]. "
                        f"Please update the model configuration or tuning bounds."
                    )
    
    def _objective_wrapper(
        self,
        params: np.ndarray,
        objective_func: Callable,
        param_names: List[str],
        evaluation_callback: Optional[Callable] = None
    ) -> float:
        """
        Wrapper for objective function that tracks history and enforces bounds.
        
        Args:
            params: Parameter values
            objective_func: Objective function to call
            param_names: Names of parameters
            evaluation_callback: Optional callback function for each function evaluation
            
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
            'evaluation': len(self.history),
            'objective': obj_value,
            'parameters': dict(zip(param_names, params.tolist()))  # Convert to list for JSON serialization
        }
        self.history.append(history_entry)
        
        # Update best
        if self.best_value is None or obj_value < self.best_value:
            self.best_value = obj_value
            self.best_params = params.copy()
        
        # Call evaluation callback if provided
        if evaluation_callback:
            evaluation_callback(history_entry)
        
        return obj_value
    
    def _callback_differential_evolution(
        self,
        intermediate_result: OptimizeResult,
        param_names: List[str],
        evaluation_callback: Optional[Callable]
    ) -> None:
        """
        Callback for differential_evolution when running in parallel.
        Runs in main process after each generation. Updates history with 'generation'
        key and calls evaluation_callback to print progress.
        
        Args:
            intermediate_result: OptimizeResult with .x (best params) and .fun (best objective)
            param_names: Parameter names for history entry
            evaluation_callback: Callback to invoke with history_entry (e.g. for printing)
        """
        x = np.asarray(intermediate_result.x)
        fun = float(intermediate_result.fun)
        history_entry = {
            'generation': len(self.history),
            'objective': fun,
            'parameters': dict(zip(param_names, x.tolist()))
        }
        self.history.append(history_entry)
        if self.best_value is None or fun < self.best_value:
            self.best_value = fun
            self.best_params = x.copy()
        if evaluation_callback:
            evaluation_callback(history_entry)
    
    def optimize(
        self,
        objective_func: Callable,
        param_names: List[str],
        bounds: List[Tuple[float, float]],
        x0: Optional[np.ndarray] = None,
        evaluation_callback: Optional[Callable] = None
    ) -> OptimizeResult:
        """
        Run optimization.
        
        Args:
            objective_func: Objective function that takes parameter array and returns scalar
            param_names: List of parameter names
            bounds: List of (min, max) tuples for each parameter
            x0: Initial guess (optional, required for some algorithms)
            evaluation_callback: Optional callback function for each function evaluation
            
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
        self._validate_optimization_inputs(bounds, x0, param_names)
        
        # For DE with workers>1, our callback runs in main process; skip per-eval callback in wrapper
        opts = dict(self.options)
        use_de_callback = self.algorithm == "differential_evolution" and opts.get('workers', 1) not in (0, 1)
        wrapper_eval_cb = None if use_de_callback else evaluation_callback
        
        # Create wrapped objective
        wrapped_obj = partial(
            self._objective_wrapper,
            objective_func=objective_func,
            param_names=param_names,
            evaluation_callback=wrapper_eval_cb
        )
        
        # Convert bounds to scipy format
        scipy_bounds = bounds
        
        if self.algorithm == "differential_evolution":
            de_opts = _coerce_numeric_options(self.options)
            if use_de_callback:
                user_cb = de_opts.pop('callback', None)
                def combined_cb(intermediate_result):
                    self._callback_differential_evolution(intermediate_result, param_names, evaluation_callback)
                    if user_cb:
                        user_cb(intermediate_result)
                de_opts['callback'] = combined_cb
            result = differential_evolution(
                wrapped_obj,
                bounds=scipy_bounds,
                **de_opts
            )
        
        elif self.algorithm == "Nelder-Mead":
            if x0 is None:
                # Use center of bounds as initial guess
                x0 = np.array([(b[0] + b[1]) / 2 for b in bounds])
            nm_opts = _coerce_numeric_options(self.options)
            result = minimize(
                wrapped_obj,
                x0=x0,
                method='Nelder-Mead',
                options=nm_opts
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
        """Print supported algorithms."""
        print("Supported algorithms:", ", ".join(sorted(cls.SUPPORTED_ALGORITHMS)))
        print("Use scipy's exact parameter names in YAML; invalid options will raise from scipy.")
