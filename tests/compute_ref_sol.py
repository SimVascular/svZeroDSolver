import os
import json
import pysvzerod

this_file_dir = os.path.abspath(os.path.dirname(__file__))

def compute_ref_sol(testname):
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
    with open(result_filename, 'w') as f:
        f.write(result.to_json())

    # print for confirmation
    print(f'Reference solution for test case {testname} computed and saved to {result_filename}. Please verify that the results are as expected.')
    
def compute_all_ref_sol():
    """
    Compute reference solutions for all test cases in the 'cases' directory.
    This function iterates through all JSON files in the 'cases' directory
    and computes their reference solutions.
    """
    testfiles = [f for f in os.listdir(os.path.join(this_file_dir, 'cases')) if f.endswith('.json')]
    testfiles.remove("steadyFlow_calibration.json")
    for testfile in testfiles:
        compute_ref_sol(testfile)


if __name__ == "__main__":
    # compute the reference solution for a specified test case
    import sys
    if len(sys.argv) != 2:
        print("Usage: python compute_ref_sol.py <test_case_name>")
        sys.exit(1)
    else:
        if sys.argv[1] == "all":
            compute_all_ref_sol()
        else:
            # compute reference solution for a single test case
            print('computing reference solution for test case:', sys.argv[1])
            test_case_name = sys.argv[1]
            compute_ref_sol(test_case_name)