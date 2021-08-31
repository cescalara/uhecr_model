#!/bin/bash

# directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Initialization of arrays
sources=('2FHL_250Mpc')
# sources=('SBG_23')
detectors=('TA2015')
models=('joint' 'joint_gmf')
# sim_models=('joint' 'joint_gmf')
# dtypes=('sim' 'real')
seeds=('19990308' '16852056')
ptypes=('p' 'N')
# ptypes=('Si' 'Fe')

# for model in "${models[@]}"; 
# do
#     for ptype in "${ptypes[@]}";
#     do
#         for seed in "${seeds[@]}";
#         do
#             printf "Config :%s\n\n" "$src joint_gmf data $ptype $seed"
#             python $SCRIPT_DIR/fit_model.py --source $src --detector $detector --model $model --dtype real --seed $seed --ptype $ptype --tight_B
#         done
#     done
# done

# do all composition-related fitting
# this simply means to run with tightening priors

# gmf fits to data, only using joint+gmf model
# iterate over all ptypes
for src in "${sources[@]}";
do
    for detector in "${detectors[@]}";
    do
        for ptype in "${ptypes[@]}"; 
        do
            for model in "${models[@]}"; 
            do
                for seed in "${seeds[@]}";
                do
                    printf "Config :%s\n\n" "$src $model data $ptype $seed"
                    python $SCRIPT_DIR/fit_model.py --source $src --detector $detector --model $model --dtype real --seed $seed --ptype $ptype
                    python $SCRIPT_DIR/fit_model.py --source $src --detector $detector --model $model --dtype real --seed $seed --ptype $ptype --tight_B
                done
            done
        done
    done
done
