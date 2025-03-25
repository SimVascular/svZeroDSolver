import os
import json
import pysvzerod
import argparse

def compute_ref_sol(testname):
    '''
    compute reference solution for a test case

    :param testname: name of the test case json to compute reference solution for
    '''
    # compute result
    result = pysvzerod.simulate(json.load(open(os.path.join('cases', testname))))

    # save result
    result_filename = os.path.join('cases', 'results', 'result_' + testname)

    # save to json
    with open(result_filename, 'w') as f:
        f.write(result.to_json())

    # print for confirmation
    print(f'Reference solution for test case {testname} computed and saved to {result_filename}. Please verify that the results are as expected.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute reference solution for a test case')
    parser.add_argument('testname', help='name of the test case json file')
    args = parser.parse_args()
    compute_ref_sol(args.testname)