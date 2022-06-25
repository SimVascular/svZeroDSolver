import argparse
from svzerodsolver.__init__ import run_simulation_from_config
import json
import os
import numpy as np


def main():

    parser = argparse.ArgumentParser(description="svZeroDSolver")

    parser.add_argument("input_file", help="Path to 0d solver input file.")
    parser.add_argument("output_file", help="Path to 0d solver output file.")

    parser.add_argument(
        "-sic",
        "--steady_ic",
        type=bool,
        default=True,
        help=(
            "Run the pulsatile 0d simulation using the steady-state solution from "
            "the equivalent steady 0d model as the initial conditions."
        ),
    )

    args = parser.parse_args()

    with open(args.input_file) as input_file:
        config = json.load(input_file)

    zero_d_results_branch = run_simulation_from_config(
        parameters=config, use_steady_soltns_as_ics=args.steady_ics
    )

    with open(args.output_file, "w") as ff:
        json.dump(zero_d_results_branch, ff, indent=4)
