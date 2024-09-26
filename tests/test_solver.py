import os
import json
import numpy as np
import os
import json
import pandas as pd
import pytest

import sys
sys.path.append(os.path.dirname(__file__))

from .utils import run_test_case_by_name, get_result, run_with_reference, RTOL_FLOW, RTOL_PRES

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
                                      'pulsatileFlow_R_RCR_mean_derivative_variable.json'
                                      ])
def test_solver(testfile):
    '''
    run all test cases and compare against stored reference solution
    '''

    this_file_dir = os.path.abspath(os.path.dirname(__file__))

    results_dir = os.path.join(this_file_dir, 'cases', 'results')

    ref = pd.read_json(os.path.join(results_dir, f'result_{testfile}'))

    run_with_reference(ref, os.path.join(this_file_dir, 'cases', testfile))
