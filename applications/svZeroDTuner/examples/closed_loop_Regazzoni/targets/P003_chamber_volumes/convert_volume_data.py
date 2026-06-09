#!/usr/bin/env python3
"""
Convert chamber volume CSVs (with RR% and mL) into time–series target files
(time in seconds, volume in m^3) suitable for svZeroDTuner time_series targets.

INPUT FILES:
- volume.csv          with columns: RR%, V_LA, V_LV, V_RA, V_RV  (volumes in mL)
- la_volume_manual.csv with columns: RR%, V_LA                  (volume in mL)

OUTPUT FILES (in the same directory as this script):
- target_V_LA.csv
- target_V_LV.csv
- target_V_RA.csv
- target_V_RV.csv
- target_V_LA_manual.csv

Each output file has columns: time (seconds), value (m^3)
"""

import os
import pandas as pd
from typing import Dict

# ============================================================
# CONFIG: EDIT THESE VALUES AS NEEDED
# ============================================================

# Cardiac period (duration of one RR interval) in seconds
CARDIAC_PERIOD_SEC = 0.689

# Values from ECG for a particular patient
# PR interval in seconds
PR_INTERVAL_SEC = 0.182

# QRS duration in seconds
QRS_DURATION_SEC = 0.088

# Input filenames (relative to this script's directory)
VOLUME_FILE = "volume.csv"
LA_MANUAL_FILE = "la_volume_manual.csv"

# Mapping from input column name -> output CSV filename
VOLUME_COLUMNS_MODEL: Dict[str, str] = {
    "V_LV": "target_V_LV.csv",
    "V_RA": "target_V_RA.csv",
    "V_RV": "target_V_RV.csv",
}

VOLUME_COLUMNS_MANUAL: Dict[str, str] = {
    "V_LA": "target_V_LA.csv",
}

# ============================================================
# HELPER FUNCTIONS
# ============================================================


def rr_percent_to_time(rr_percent_series, cardiac_period_sec: float):
    """
    Convert RR% (0–100 over one cardiac cycle) to time in seconds.
    RR%=0 corresponds to the R-wave. Assuming the start of the P-wave corresponds
    to t=0, then RR%=0 corresponds to the PR interval + the QRS interval/2.
    time = (rr/100) * cardiac_period + PR interval + QRS interval/2, wrapped into [0, cardiac_period).
    """
    rr = rr_percent_series.astype(float) / 100.0
    time = rr * cardiac_period_sec + PR_INTERVAL_SEC + QRS_DURATION_SEC / 2
    return time % cardiac_period_sec


def convert_volume_file(
    csv_path: str,
    column_to_output: Dict[str, str],
):
    """
    Convert a volume CSV (RR% and mL) into time–series target CSVs (time [s], value [m^3]).
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Input file not found: {csv_path}")

    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()  # handle " V_LA" -> "V_LA"

    if "RR%" not in df.columns:
        raise ValueError(f"'RR%' column not found in {csv_path}. Columns: {list(df.columns)}")

    times = rr_percent_to_time(df["RR%"], CARDIAC_PERIOD_SEC)
    directory = os.path.dirname(os.path.abspath(csv_path))

    for col, out_name in column_to_output.items():
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in {csv_path}. Columns: {list(df.columns)}")

        values_m3 = df[col].astype(float) * 1e-6  # mL -> m^3

        out_df = (
            pd.DataFrame({"time": times, "value": values_m3})
            .sort_values("time")
        )

        out_path = os.path.join(directory, out_name)
        out_df.to_csv(out_path, index=False)
        print(f"Saved {out_path} ({len(out_df)} points)")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    convert_volume_file(
        csv_path=os.path.join(script_dir, VOLUME_FILE),
        column_to_output=VOLUME_COLUMNS_MODEL
    )

    convert_volume_file(
        csv_path=os.path.join(script_dir, LA_MANUAL_FILE),
        column_to_output=VOLUME_COLUMNS_MANUAL
    )


if __name__ == "__main__":
    main()
