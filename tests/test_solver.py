import os
import json
import numpy as np
import os
import json
import pandas as pd
import pytest

import sys
sys.path.append(os.path.dirname(__file__))

from .utils import run_with_reference, RTOL_PRES, RTOL_FLOW

EXPECTED_FAILURES = {
    'closedLoopHeart_singleVessel_mistmatchPeriod.json',
    'pulsatileFlow_R_RCR_mismatchPeriod.json'
}

@pytest.mark.parametrize("testfile", ['chamber_elastance_inductor.json', 
                                      'steadyFlow_R_R.json', 
                                      'pulsatileFlow_R_coronary_cycle_error.json', 
                                      'pulsatileFlow_R_coronary.json', 
                                      'coupledBlock_closedLoopHeart_withCoronaries.json', 
                                      'pulsatileFlow_R_RCR_derivative.json', 
                                      'steadyFlow_bifurcationR_R1.json', 
                                      'closedLoopHeart_singleVessel.json', 
                                      'pulsatileFlow_R_RCR_variable.json', 
                                      'coupledBlock_closedLoopHeart_singleVessel.json', 
                                      'steadyFlow_R_steadyPressure.json', 
                                      'steadyFlow_bifurcationR_R1_blockNames.json', 
                                      'pulsatileFlow_R_RCR.json', 
                                      'closedLoopHeart_withCoronaries.json', 
                                      'steadyFlow_RL_R.json', 
                                      'steadyFlow_bifurcationR_R2.json', 
                                      'pulsatileFlow_R_RCR_derivative_variable.json', 
                                      'steadyFlow_confluenceR_R.json', 
                                      'steadyFlow_R_RCR.json', 
                                      'steadyFlow_stenosis_R.json',
                                      'pulsatileFlow_R_RCR_mean.json', 
                                      'pulsatileFlow_R_RCR_mean_derivative.json', 
                                      'pulsatileFlow_CStenosis_steadyPressure.json', 
                                      'pulsatileFlow_R_RCR_mean_variable.json', 
                                      'steadyFlow_RC_R.json', 
                                      'steadyFlow_R_coronary.json', 
                                      'steadyFlow_RLC_R.json', 
                                      'steadyFlow_blood_vessel_junction.json', 
                                      'valve_tanh.json', 
                                      'pulsatileFlow_bifurcationR_RCR_cycle_error.json', 
                                      'pulsatileFlow_R_RCR_mean_derivative_variable.json',
                                      'closedLoopHeart_singleVessel_mistmatchPeriod.json',
                                      'pulsatileFlow_R_RCR_mismatchPeriod.json',
                                      'pulsatileFlow_CStenosis_steadyPressure_definedPeriod.json',
                                      'chamber_sphere.json',
                                      'piecewise_Chamber_and_Valve.json'
                                      ])
def test_solver(testfile):
    '''
    run all test cases and compare against stored reference solution
    '''

    # default tolerances
    rtol_pres = RTOL_PRES
    rtol_flow = RTOL_FLOW
    if 'coupledBlock_closedLoopHeart_withCoronaries.json' in testfile:
        rtol_pres = 2.0e-1
        rtol_flow = 2.0e-1

    this_file_dir = os.path.abspath(os.path.dirname(__file__))

    results_dir = os.path.join(this_file_dir, 'cases', 'results')

    if testfile in EXPECTED_FAILURES:
        pytest.xfail(reason=f"Known failure for test case: {testfile}")

    ref = pd.read_json(os.path.join(results_dir, f'result_{testfile}'))

    run_with_reference(ref, os.path.join(this_file_dir, 'cases', testfile), rtol_pres, rtol_flow)
