"""
Output Extractor for sv0D Tuning Framework

Extracts outputs from sv0D simulation results using exact variable names from results.csv.
Supports time_series, min, max, and mean extraction types.
"""

import numpy as np
import pandas as pd
from typing import Union, Literal, Optional
import pysvzerod


class OutputExtractor:
    """
    Extracts outputs from sv0D simulation results.
    
    Uses exact variable names from sv0D results (e.g., "flow:AR_SYS:J0").
    Supports extraction of time series or scalar metrics (min, max, mean).
    """
    
    def __init__(self, solver: pysvzerod.Solver):
        """
        Initialize output extractor with a solver instance.
        
        Args:
            solver: pysvzerod.Solver instance (must have run() called)
        """
        self.solver = solver
        self._result_df = None
        self._times = None
    
    def _ensure_results(self):
        """Ensure simulation has been run and results are available."""
        if self._result_df is None:
            self._result_df = self.solver.get_full_result()
            self._times = self.solver.get_times()
    
    def extract(
        self, 
        output_name: str, 
        extraction_type: Literal["time_series", "min", "max", "mean"] = "time_series"
    ) -> Union[np.ndarray, float]:
        """
        Extract output from simulation results.
        
        Args:
            output_name: Exact output name from sv0D results (e.g., "flow:AR_SYS:J0")
            extraction_type: Type of extraction:
                - "time_series": Return full time series array
                - "min": Return minimum value
                - "max": Return maximum value
                - "mean": Return mean value
                
        Returns:
            Time series array (if extraction_type="time_series") or scalar value
        """
        self._ensure_results()
        
        # Get the time series for this output
        try:
            time_series = self.solver.get_single_result(output_name)
        except Exception as e:
            raise ValueError(f"Output '{output_name}' not found in simulation results: {e}")
        
        # Apply extraction type
        if extraction_type == "time_series":
            return np.array(time_series)
        elif extraction_type == "min":
            return float(np.min(time_series))
        elif extraction_type == "max":
            return float(np.max(time_series))
        elif extraction_type == "mean":
            return float(np.mean(time_series))
        else:
            raise ValueError(f"Unknown extraction_type: {extraction_type}")
    
    def get_times(self) -> np.ndarray:
        """
        Get time array from simulation.
        
        Returns:
            Time array
        """
        self._ensure_results()
        return np.array(self._times)
    
    def get_all_output_names(self) -> list:
        """
        Get list of all available output names.
        
        Returns:
            List of output names
        """
        self._ensure_results()
        return list(self._result_df['name'].unique())
    
    def extract_multiple(
        self, 
        outputs: list[dict]
    ) -> dict:
        """
        Extract multiple outputs at once.
        
        Args:
            outputs: List of dicts with keys:
                - name: Output name
                - type: Extraction type (time_series, min, max, mean)
                
        Returns:
            Dictionary mapping output names to extracted values
        """
        results = {}
        for output_spec in outputs:
            name = output_spec['name']
            extraction_type = output_spec.get('type', 'time_series')
            results[name] = self.extract(name, extraction_type)
        return results
