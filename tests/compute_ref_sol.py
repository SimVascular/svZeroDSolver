import os
import json
import pysvzerod

def compute_ref_sol(testname):
    '''
    compute reference solution for a test case

    :param testname: name of the test case json to compute reference solution for
    '''

    # testfiles = [f for f in os.listdir('tests/cases') if os.path.isfile(os.path.join('tests/cases', f))]

    # we only want to test the solver, not the calibrator
    # testfiles.remove("steadyFlow_calibration.json")

    # compute result
    result = pysvzerod.simulate(json.load(open(os.path.join('cases', testname))))

    # save result
    result_filename = os.path.join('cases/results', 'result_' + testname)

    # save to json
    with open(result_filename, 'w') as f:
        f.write(result.to_json())

    # print for confirmation
    print(f'Reference solution for test case {testname} computed and saved to {result_filename}. Please verify that the results are as expected.')
    
if __name__ == "__main__":
    # compute the reference solution for a specified test case
    import sys
    if len(sys.argv) != 2:
        print("Usage: python compute_ref_sol.py <test_case_name>")
        sys.exit(1)
    else:
        print('computing reference solution for test case:', sys.argv[1])
        test_case_name = sys.argv[1]
        compute_ref_sol(test_case_name)