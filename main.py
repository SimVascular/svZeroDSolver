import argparse
from svzerodsolver.__init__ import run_simulation_from_config
import json
import os
import numpy as np

parser = argparse.ArgumentParser(description="svZeroDSolver")

parser.add_argument("input_file", help="Path to 0d solver input file.")

parser.add_argument(
    "-sic",
    "--steady_ic",
    action="store_true",
    help=(
        "Run the pulsatile 0d simulation using the steady-state solution from "
        "the equivalent steady 0d model as the initial conditions."
    ),
)

args = parser.parse_args()

with open(args.input_file) as input_file:
    config = json.load(input_file)

zero_d_results_branch, zero_d_results_all = run_simulation_from_config(
    parameters=config,
    use_steady_soltns_as_ics=args.steady_ic,
)

# postprocessing
zero_d_input_file_name = os.path.splitext(args.input_file)[0]
if args.saveAll:
    np.save(zero_d_input_file_name + "_all_results", zero_d_results_all)

if args.saveBranch:
    np.save(zero_d_input_file_name + "_branch_results", zero_d_results_branch)
