"""
Shared expression evaluation for targets and QoIs.
Expressions reference simulation output names (e.g. Vc:LV, pressure:AV:AR_SYS)
and are evaluated with numpy (np) and output_arrays["output_name"] -> array.
"""

import numpy as np
from typing import Callable, Dict, List, Literal, Tuple, Union


class Expression:
    """
    Represents an expression that references simulation outputs (e.g. np.max(Vc:LV)).
    Each target or QoI has an Expression object with an evaluate method.

    Security note: expressions are compiled and executed via exec() with full Python
    builtins available. Config files should be treated as trusted code — do not run
    configs from untrusted sources.
    """

    def __init__(
        self,
        expression_str: str,
        kind: Literal["scalar", "time_series"] = "scalar",
    ):
        """
        Args:
            expression_str: The expression string (e.g. np.max(Vc:LV)).
            kind: "scalar" returns float; "time_series" returns (values_array, times_array).
        """
        self.expression_str = expression_str
        self.kind = kind
        # Cache: frozenset(available_outputs) -> (output_names, eval_func)
        self._compile_cache: Dict[frozenset, Tuple[List[str], Callable]] = {}

    def __getstate__(self):
        """
        Make Expression pickle-safe for multiprocessing by dropping compiled callables.
        They are recreated lazily on first evaluate() in each worker process.
        """
        state = self.__dict__.copy()
        state["_compile_cache"] = {}
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        if "_compile_cache" not in self.__dict__:
            self._compile_cache = {}

    def output_names(self, available_outputs: List[str]) -> List[str]:
        """Return output names referenced in this expression (for extracting data)."""
        return [
            out_name
            for out_name in available_outputs
            if out_name in self.expression_str
        ]

    def _get_compiled(
        self, available_outputs: List[str]
    ) -> Tuple[List[str], Callable[[Dict[str, np.ndarray]], object]]:
        """
        Return (output_names, eval_func) for this expression and available_outputs.
        Eval func takes output_arrays and returns the expression result.
        Compilation is cached per set of available_outputs.
        """
        key = frozenset(available_outputs)
        if key not in self._compile_cache:
            output_names = [
                o for o in available_outputs if o in self.expression_str
            ]
            eval_str = self.expression_str
            for out_name in sorted(output_names, key=len, reverse=True):
                eval_str = eval_str.replace(
                    out_name, f'output_arrays["{out_name}"]'
                )
            namespace: Dict[str, object] = {"np": np}
            code = compile(
                f"def _eval(output_arrays):\n    return {eval_str}",
                "<expr>",
                "exec",
            )
            exec(code, namespace)
            self._compile_cache[key] = (output_names, namespace["_eval"])  # type: ignore[index]
        return self._compile_cache[key]

    def evaluate(
        self,
        unique_outputs: Dict[str, Dict],
        available_outputs: List[str],
    ) -> Union[float, Tuple[np.ndarray, np.ndarray]]:
        """
        Evaluate the expression with the given outputs.
        Returns float for scalar, (values_array, times_array) for time_series.
        """
        if self.kind == "scalar":
            return self._evaluate_scalar(unique_outputs, available_outputs)
        else:
            return self._evaluate_time_series(unique_outputs, available_outputs)

    def _build_output_arrays(
        self,
        unique_outputs: Dict[str, Dict],
        output_names: List[str],
    ) -> Dict[str, np.ndarray]:
        """Build dict mapping output name -> array for the given output names."""
        output_arrays: Dict[str, np.ndarray] = {}
        for out_name in output_names:
            if out_name not in unique_outputs:
                continue
            data = unique_outputs[out_name]
            arr = data.get("time_series", data.get("value"))
            if isinstance(arr, (int, float)):
                arr = np.array([arr])
            output_arrays[out_name] = np.asarray(arr)
        return output_arrays

    def _evaluate_scalar(
        self,
        unique_outputs: Dict[str, Dict],
        available_outputs: List[str],
    ) -> float:
        output_names, eval_func = self._get_compiled(available_outputs)
        output_arrays = self._build_output_arrays(unique_outputs, output_names)
        try:
            result = eval_func(output_arrays)
        except Exception as e:
            raise ValueError(f"Expression evaluation failed: {e}") from e
        arr = np.asarray(result)
        if arr.size != 1:
            raise ValueError(
                f"Scalar expression must yield a single value, "
                f"got shape {arr.shape} (size {arr.size})"
            )
        return float(arr.item())

    def _evaluate_time_series(
        self,
        unique_outputs: Dict[str, Dict],
        available_outputs: List[str],
    ) -> Tuple[np.ndarray, np.ndarray]:
        output_names, eval_func = self._get_compiled(available_outputs)
        output_arrays = self._build_output_arrays(unique_outputs, output_names)
        try:
            result = eval_func(output_arrays)
        except Exception as e:
            raise ValueError(f"Expression evaluation failed: {e}") from e
        arr = np.asarray(result)
        if arr.ndim != 1:
            raise ValueError(
                f"Time series expression must return 1D array, got shape {arr.shape}"
            )
        times = None
        for out_name in output_arrays:
            if out_name in unique_outputs and "times" in unique_outputs[out_name]:
                times = np.asarray(unique_outputs[out_name]["times"])
                break
        if times is None or len(times) != len(arr):
            times = np.arange(len(arr))
        return arr, times
