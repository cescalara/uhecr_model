#!/bin/bash

# directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Initialization of arrays
sources=('SBG_23' '2FHL_250Mpc' 'swift_BAT_213')
detectors=('TA2015' 'auger2014')
models=('arrival_direction' 'joint')
sim_models=('joint' 'joint_gmf')
dtypes=('sim' 'real')
seeds=('19990308' '16852056' '65492186' '09953497' '9999999')
ptypes=('p' 'He' 'N' 'O' 'Si' 'Fe')

# do all gmf-related sim + fitting

# gmf sims
# for ptype in "${ptypes[@]}";
# do
#     for src in "${sources[@]}"; 
#     do
#         for seed in "${seeds[@]}"; 
#         do
#             printf  "Config :%s\n\n" "$ptype $src $seed"
#             python $SCRIPT_DIR/simulate_data.py --source $src --detector TA2015 --model joint_gmf --seed $seed --ptype $ptype
#         done
#     done
# done


# gmf fits to joint simulation
# models only use arrival and joint

# for src in "${sources[@]}";
#     do
#         for model in "${models[@]}"; 
#         do
#             for seed in "${seeds[@]}";
#             do
#                 printf "Config :%s\n\n" "$src sim $model $seed"
#                 python $SCRIPT_DIR/fit_model.py --source $src --detector TA2015 --model $model --sim_model joint --dtype sim --seed $seed --ptype p
#             done
#         done
#     done
# done

# gmf fits to data
# use only arrival and joint, since ptype will always only be p
for src in "${sources[@]}";
    do
        for model in "${models[@]}"; 
        do
            for seed in "${seeds[@]}";
            do
                printf "Config :%s\n\n" "$src data $model $seed"
                python $SCRIPT_DIR/fit_model.py --source $src --detector TA2015 --model $model --dtype real --seed $seed --ptype p
            done
        done
    done
done