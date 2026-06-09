"""
Optimizer wrappers for sv0D Tuning Framework

Supports scipy optimizers and parallelization.
"""

import numpy as np
from typing import Callable, List, Tuple, Dict, Optional, Any
from scipy.optimize import minimize, differential_evolution, OptimizeResult


# Small tolerance for "objective hits zero" (floating-point safety)
_OBJECTIVE_ZERO_TOL = 1e-12

# Tolerance for "close to bound" warning: within this fraction of the range from a bound
_BOUND_TOLERANCE = 0.01


class _ScaledObjective:
    """Pickle-safe objective wrapper for scaled optimizer-space parameters."""

    def __init__(
        self,
        objective_func: Callable[[np.ndarray], float],
        to_phys: Callable[[np.ndarray], np.ndarray],
    ):
        self.objective_func = objective_func
        self.to_phys = to_phys

    def __call__(self, x_opt: np.ndarray) -> float:
        return self.objective_func(self.to_phys(np.asarray(x_opt)))


def _check_params_near_bounds(
    params: Dict[str, float], parameters: List[Dict]
) -> List[tuple]:
    """Return (name, value, bound_type, bound_value) for params near bounds."""
    near_bounds = []
    for p in parameters:
        name = p['name']
        if name not in params:
            continue
        bounds = p.get('bounds')
        if not bounds or len(bounds) != 2:
            continue
        lower, upper = float(bounds[0]), float(bounds[1])
        value = float(params[name])
        range_val = upper - lower
        if range_val <= 0:
            continue
        tol = _BOUND_TOLERANCE * range_val
        if value <= lower + tol:
            near_bounds.append((name, value, 'min', lower))
        elif value >= upper - tol:
            near_bounds.append((name, value, 'max', upper))
    return near_bounds


class ObjectiveReachedZero(Exception):
    """Raised when optimization terminates because objective reached zero (Nelder-Mead)."""
    pass


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
        self.options = dict(options)
        self.terminate_at_zero = self.options.pop('terminate_at_zero', True)
        self._validate_config()
        self.history = []
        self.best_value = None
        self.best_params = None
    
    def _validate_config(self):
        if self.algorithm not in self.SUPPORTED_ALGORITHMS:
            raise ValueError(
                f"Algorithm '{self.algorithm}' is not supported. "
                f"Supported: {', '.join(sorted(self.SUPPORTED_ALGORITHMS))}"
            )
    
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
    
    def master_callback(
        self,
        intermediate_result: OptimizeResult,
        param_names: List[str],
        parameters: Optional[List[Dict]] = None,
    ) -> Optional[bool]:
        """
        Callback function for optimization algorithms. Update history, best, print progress, warn near bounds.
        Checks if objective reached zero and stops optimization if necessary (returns True for DE, or raises ObjectiveReachedZero for NM, otherwise None).
        """
        # Extract parameters and objective value from OptimizeResult
        x = np.asarray(intermediate_result.x)
        fun = float(intermediate_result.fun)

        # Convert x to physical space if scaling is used
        if self._use_scaling:
            x = self._scale_to_phys(x)
           
        params_dict = dict(zip(param_names, x.tolist()))
        
        # Update history and best parameters
        history_entry = {
            'iteration': len(self.history),
            'objective': fun,
            'parameters': params_dict,
        }
        self.history.append(history_entry)
        if self.best_value is None or fun < self.best_value:
            self.best_value = fun
            self.best_params = x.copy()

        # Print progress
        step = len(self.history) - 1
        print(f"Iteration {step:3d}: Objective = {fun:.6e}", end="")
        if len(params_dict) <= 3:
            param_str = ", ".join([f"{name}={val:.3e}" for name, val in params_dict.items()])
            print(f" | {param_str}")
        else:
            print()

        # Warn if any parameter is close to its bounds
        if parameters:
            for name, value, bound_type, bound_value in _check_params_near_bounds(params_dict, parameters):
                print(f"  WARNING: {name}={value:.3e} is near {bound_type} bound ({bound_value:.3e})")

        # If terminating optimization at zero objective, return True for DE, raise ObjectiveReachedZero for NM, otherwise None.
        if self.terminate_at_zero and fun <= _OBJECTIVE_ZERO_TOL:
            if self.algorithm == "differential_evolution":
                return True
            raise ObjectiveReachedZero("Objective reached zero; stopping optimization")
        return None

    def optimize(
        self,
        objective_func: Callable,
        param_names: List[str],
        bounds: List[Tuple[float, float]],
        x0: Optional[np.ndarray] = None,
        parameters: Optional[List[Dict]] = None,
        param_scaling_to_opt_space: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        param_scaling_to_phys_space: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    ) -> OptimizeResult:
        """
        Run optimization.

        Bounds and x0 are in physical space. If scaling callables are
        provided, the optimizer works in scaled/optimizer space and returns result.x and
        get_best() in physical space.

        Args:
            objective_func: Objective function that takes physical-space parameter
                array and returns scalar.
            param_names: List of parameter names.
            bounds: List of (min, max) tuples in physical space.
            x0: Initial guess in physical space (optional).
            parameters: Optional list of param dicts with 'name' and 'bounds'
                (used for near-bounds warnings).
            param_scaling_to_opt_space: Optional; convert physical -> scaled/optimizer space.
            param_scaling_to_phys_space: Optional; convert scaled/optimizer -> physical space.

        Returns:
            OptimizeResult with result.x in physical space.
        """
        self.history = []
        self.best_value = None
        self.best_params = None
        self._validate_optimization_inputs(bounds, x0, param_names)

        # Ensure numeric optimization options are numeric types.
        opts = _coerce_numeric_options(self.options)

        # Callback function supporting both SciPy styles:
        # - callback(intermediate_result=OptimizeResult)
        # - callback(x, convergence=...)
        def callback(*args, **kwargs):
            if "intermediate_result" in kwargs:
                return self.master_callback(
                    kwargs["intermediate_result"], param_names, parameters
                )

            if len(args) == 1 and isinstance(args[0], OptimizeResult):
                return self.master_callback(args[0], param_names, parameters)

            x = None
            if len(args) >= 1:
                x = np.asarray(args[0])
            elif "xk" in kwargs:
                x = np.asarray(kwargs["xk"])
            elif "x" in kwargs:
                x = np.asarray(kwargs["x"])

            if x is None:
                return None

            # Legacy DE callback doesn't provide objective value; evaluate it here.
            x_phys = self._scale_to_phys(x) if self._use_scaling else x
            fun = float(objective_func(x_phys))
            return self.master_callback(
                OptimizeResult(x=x, fun=fun), param_names, parameters
            )

        # Use center of bounds as initial guess if no initial guess is provided
        if x0 is None:
            x0 = np.array([(b[0] + b[1]) / 2 for b in bounds])
        
        # Handle scaling of parameter values
        self._use_scaling = (
            param_scaling_to_opt_space is not None
            and param_scaling_to_phys_space is not None
        )
        if self._use_scaling:
            self._scale_to_opt = param_scaling_to_opt_space
            self._scale_to_phys = param_scaling_to_phys_space

            # Scale bounds
            bounds_arr = np.array(bounds)
            bounds = list(zip(
                self._scale_to_opt(bounds_arr[:, 0]),
                self._scale_to_opt(bounds_arr[:, 1]),
            ))

            # Scale initial guess
            x0 = self._scale_to_opt(np.asarray(x0))

            # Objective function in optimizer space (module-level callable for multiprocessing compatibility)
            obj_fun = _ScaledObjective(objective_func, self._scale_to_phys)
        else:
            obj_fun = objective_func

        if self.algorithm == "differential_evolution":
            opts['callback'] = callback
            result = differential_evolution(obj_fun, bounds=bounds, **opts)
        elif self.algorithm == "Nelder-Mead":
            try:
                result = minimize(
                    obj_fun, x0=x0, method='Nelder-Mead',
                    bounds=bounds, options=opts, callback=callback
                )
            except ObjectiveReachedZero:
                # Terminated early because objective reached zero
                result = OptimizeResult(
                    # best_params is in physical space from callback, so convert to optimizer space
                    # to be consistent with the optimizer's returned result
                    x=self._scale_to_opt(self.best_params) if self._use_scaling else self.best_params,
                    success=True,
                    fun=self.best_value,
                    message="Optimization terminated: objective reached zero",
                    nfev=len(self.history),
                )

        # Use the optimizer's returned result for best value/params (not best seen during optimization)
        if hasattr(result, 'x') and hasattr(result, 'fun'):
            if self._use_scaling:
                result.x = self._scale_to_phys(np.asarray(result.x))
            self.best_params = result.x
            self.best_value = result.fun

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
