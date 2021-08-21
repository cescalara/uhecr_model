#!/bin/bash
# 
# Initialization of arrays
sources = ('SBG_23', '2FHL_250Mpc', 'swift_BAT_213')
detectors = ('TA2015', 'auger2014')
models = ('arrival_direction', 'joint', 'joint_gmf')
dtype = ('sim', 'real')
seeds = ('19990308', '16852056', '09953497', '9999999')
ptypes = ('p', 'He', 'Li', 'C', 'N', 'O')
#
python ./simulate_data.py 