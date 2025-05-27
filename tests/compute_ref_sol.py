import os
import json
import pysvzerod
import argparse
import sys

this_file_dir = os.path.abspath(os.path.dirname(__file__))

def compute_reference_solution(testname):
    '''
    compute reference solution for a test case

    :param testname: name of the test case json to compute reference solution for
    '''

    # testfiles = [f for f in os.listdir('tests/cases') if os.path.isfile(os.path.join('tests/cases', f))]

    # we only want to test the solver, not the calibrator
    # testfiles.remove("steadyFlow_calibration.json")

    # compute result
    result = pysvzerod.simulate(json.load(open(os.path.join(this_file_dir, 'cases', testname))))

    # save result
    result_filename = os.path.join(this_file_dir, 'cases', 'results', 'result_' + testname)

    # save to json
    result.to_json(result_filename)

    print(f'Reference solution for test case {testname} computed and saved to {result_filename}. Please verify that the results are as expected.')
    
def compute_all_ref_sol(tolerance=None):
    """
    Compute reference solutions for all test cases in the 'cases' directory.
    This function iterates through all JSON files in the 'cases' directory
    and computes their reference solutions.
    """
    testfiles = [f for f in os.listdir(os.path.join(this_file_dir, 'cases')) if f.endswith('.json')]
    testfiles.remove("steadyFlow_calibration.json")
    for testfile in testfiles:
        if tolerance is not None:
            set_absolute_tolerance(testfile, tolerance)

        compute_reference_solution(testfile)


def set_absolute_tolerance(case_name, tolerance):
    """
    Set the absolute tolerance for a specific test case.

    Args:
        case_name: Name of the test case.
        tolerance: Absolute tolerance value.
    """
    # Load the existing configuration
    with open(os.path.join(this_file_dir, "cases", case_name), "r") as f:
        config = json.load(f)

    # Set the absolute tolerance
    config["simulation_parameters"]["absolute_tolerance"] = tolerance

    # Save the updated configuration
    with open(os.path.join(this_file_dir, "cases", case_name), "w") as f:
        json.dump(config, f, indent=4)
        

def main():
    parser = argparse.ArgumentParser(description="Compute reference solution for a specified test case.")
    parser.add_argument(
        "test_case_name",
        help="Name of the test case or 'all' to compute for all test cases"
    )
    parser.add_argument(
        "--abs-tol",
        type=float,
        help="Absolute tolerance to use for the test case (only applies to single test cases)"
    )

    args = parser.parse_args()

    if args.test_case_name == "all":
        if args.abs_tol is not None:
            compute_all_ref_sol(args.abs_tol)
        else:
            compute_all_ref_sol()
    else:
        print("Computing reference solution for test case:", args.test_case_name)
        if args.abs_tol is not None:
            set_absolute_tolerance(args.test_case_name, args.abs_tol)
        compute_reference_solution(args.test_case_name)

if __name__ == "__main__":
    
    main()