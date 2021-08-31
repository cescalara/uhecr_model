#!/bin/bash

# directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Initialization of arrays
sources=('SBG_23' '2FHL_250Mpc' 'swift_BAT_213')
detectors=('TA2015' 'auger2014')
models=('arrival_direction' 'joint' 'joint_gmf')
sim_models=('joint' 'joint_gmf')
dtypes=('sim' 'real')
seeds=('19990308' '16852056')
# ptypes=('p' 'He' 'N' 'O' 'Si' 'Fe')
ptypes=('Si' 'Fe')

# do all gmf-related sim + fitting

# gmf sims
# for ptype in "${ptypes[@]}";
# do
#     for src in "${sources[@]}"; 
#     do
#         for seed in "${seeds[@]}"; 
#         do
#             printf  "Config :%s\n\n" "$ptype $src $seed"
#             python $SCRIPT_DIR/simulate_data.py --source $src --detector TA2015 --model joint_gmf --seed $seed --ptype $ptype --sim_inputs 0.5 1 3.0
#         done
#     done
# done


# gmf fits to joint+gmf simulation
# account for all models and all ptypes
# for ptype in "${ptypes[@]}"; 
# do
#     for src in "${sources[@]}";
#     do
#         for model in "${models[@]}"; 
#         do
#             for seed in "${seeds[@]}";
#             do
#                 printf "Config :%s\n\n" "$src $model sim $ptype $seed"
#                 python $SCRIPT_DIR/fit_model.py --source $src --detector TA2015 --model $model --sim_model joint_gmf --dtype sim --seed $seed --ptype $ptype
#             done
#         done
#     done
# done

# gmf fits to data, only using joint+gmf model
# iterate over all ptypes
for src in "${sources[@]}"; 
do
    for ptype in "${ptypes[@]}";
    do
        for seed in "${seeds[@]}";
        do
            printf "Config :%s\n\n" "$src joint_gmf data $ptype $seed"
            python $SCRIPT_DIR/fit_model.py --source $src --detector TA2015 --model joint_gmf --dtype real --seed $seed --ptype $ptype
        done
    done
done
