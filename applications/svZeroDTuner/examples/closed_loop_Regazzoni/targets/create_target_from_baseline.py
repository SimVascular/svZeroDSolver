#!/usr/bin/env python3
"""
Create a target CSV from a baseline_results.csv column.

Extracts evenly spaced samples over 1 cardiac cycle and saves as CSV.
Shows a plot of the full baseline and sampled points.

Modify the hardcoded paths and options in main() as needed.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    # Modify these paths and options as needed
    script_dir = os.path.dirname(os.path.abspath(__file__))
    baseline_path = os.path.join(script_dir, "..", "baseline_results", "baseline_results.csv")
    column = "pressure:AV:AR_SYS"
    num_samples = 20
    output_csv = os.path.join(script_dir, "target_pressure_ar_sys.csv")

    if not os.path.exists(baseline_path):
        raise FileNotFoundError(
            f"Baseline results not found: {baseline_path}\n"
            "Generate baseline results first with:\n"
            "python -c 'from main import run_baseline; run_baseline(\"model.json\")'"
        )

    df = pd.read_csv(baseline_path)
    if column not in df.columns:
        raise ValueError(
            f"Column '{column}' not found. Available: {list(df.columns)}"
        )

    times = df["time"].values
    values = df[column].values

    n = len(times)
    if n < 2:
        raise ValueError("Baseline results have too few points")
    indices = np.linspace(0, n - 1, num_samples, dtype=int)
    indices = np.unique(indices)

    time_samples = times[indices]
    value_samples = values[indices]

    pd.DataFrame({"time": time_samples, "value": value_samples}).to_csv(
        output_csv, index=False
    )
    print(f"Saved {output_csv} ({len(time_samples)} samples)")

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(times, values, "b-", alpha=0.4, label="Full baseline")
    ax.plot(
        time_samples,
        value_samples,
        "ro-",
        markersize=6,
        label=f"Target samples (n={len(time_samples)})",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(column)
    ax.set_title(f"Target: {column} ({num_samples} samples over 1 cardiac cycle)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
