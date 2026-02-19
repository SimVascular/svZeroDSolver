"""
Shared expression evaluation for targets and QoIs.
Expressions reference simulation output names (e.g. Vc:LV, pressure:AV:AR_SYS)
and are evaluated with numpy (np) and output_arrays["output_name"] -> array.
"""

import numpy as np
from typing import Dict, List, Literal, Tuple, Union


class Expression:
    """
    Represents an expression that references simulation outputs (e.g. np.max(Vc:LV)).
    Each target or QoI has an Expression object with an evaluate method.
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

    def output_names(self, available_outputs: List[str]) -> List[str]:
        """Return output names referenced in this expression (for extracting data)."""
        return [
            out_name
            for out_name in available_outputs
            if out_name in self.expression_str
        ]

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

    def _evaluate_scalar(
        self,
        unique_outputs: Dict[str, Dict],
        available_outputs: List[str],
    ) -> float:
        output_arrays = self._build_output_arrays(unique_outputs, available_outputs)
        eval_str = self._substitute_outputs(output_arrays)
        try:
            result = eval(eval_str, {"np": np, "output_arrays": output_arrays}, {})
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
        output_arrays = self._build_output_arrays(unique_outputs, available_outputs)
        eval_str = self._substitute_outputs(output_arrays)
        try:
            result = eval(eval_str, {"np": np, "output_arrays": output_arrays}, {})
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

    def _build_output_arrays(
        self,
        unique_outputs: Dict[str, Dict],
        available_outputs: List[str],
    ) -> Dict[str, np.ndarray]:
        """Build dict mapping output name -> array for names in this expression."""
        output_arrays = {}
        for out_name in self.output_names(available_outputs):
            if out_name not in unique_outputs:
                continue
            data = unique_outputs[out_name]
            arr = data.get("time_series", data.get("value"))
            if isinstance(arr, (int, float)):
                arr = np.array([arr])
            output_arrays[out_name] = np.asarray(arr)
        return output_arrays

    def _substitute_outputs(
        self, output_arrays: Dict[str, np.ndarray]
    ) -> str:
        """Replace output names in expression with output_arrays["name"]."""
        eval_str = self.expression_str
        for out_name in sorted(output_arrays.keys(), key=len, reverse=True):
            eval_str = eval_str.replace(
                out_name, f'output_arrays["{out_name}"]'
            )
        return eval_str
