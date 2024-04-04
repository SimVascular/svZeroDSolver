import os
import json
import pysvzerod

def get_all_results():
    '''
    compute results for all of the test cases and store them in the results directory
    '''

    testfiles = [f for f in os.listdir('tests/cases') if os.path.isfile(os.path.join('tests/cases', f))]

    # we only want to test the solver, not the calibrator
    testfiles.remove("steadyFlow_calibration.json")

    for file in testfiles:

        # compute result
        result = pysvzerod.simulate(json.load(open(os.path.join('tests/cases', file))))

        # save result
        with open(os.path.join('tests/cases/results', 'result_' + file), 'w') as f:
            f.write(result.to_json())


if __name__ == "__main__":
    get_all_results()
    

        
