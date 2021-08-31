#!/bin/bash

# directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Initialization of arrays
sources=('SBG_23' '2FHL_250Mpc')
detectors=('TA2015' 'auger2014')
models=('arrival_direction' 'joint' 'joint_gmf')
sim_models=('joint' 'joint_gmf')
dtypes=('sim' 'real')
seeds=('163591' '78520492' '51610230' '8946064' '9526433' '84910385' '11223344' '7498365' '196560486' '98765432' '318979')
ptypes=('Si' 'Fe')

# do all gmf-related sim + fitting

# gmf sims
# for src in "${sources[@]}"; 
# do
#     for ptype in "${ptypes[@]}";
#     do
#         for seed in "${seeds[@]}"; 
#         do
#             printf  "Config :%s\n\n" "$ptype $src $seed"
#             python $SCRIPT_DIR/simulate_data.py --source $src --detector TA2015 --model joint_gmf --seed $seed --ptype $ptype
#         done
#     done
# done

# for src in "${sources[@]}"; 
# do
#     for model in "${models[@]}"; 
#     do
#         for ptype in "${ptypes[@]}";
#         do
#             for seed in "${seeds[@]}"; 
#             do
#                 printf  "Config :%s\n\n" "$model $ptype $src $seed"
#                 python $SCRIPT_DIR/fit_model.py --source $src --detector TA2015 --model $model --dtype sim --seed $seed --ptype $ptype --sim_model joint_gmf
#                 python $SCRIPT_DIR/fit_model.py --source $src --detector TA2015 --model $model --dtype sim --seed $seed --ptype $ptype --sim_model joint_gmf --tight_B
#             done
#         done
#     done
# done



# for seed in "${seeds[@]}"; 
# do
#     printf  "Config :%s\n\n" "$seed"
#     python $SCRIPT_DIR/simulate_data.py --source SBG_23 --detector TA2015 --model joint_gmf --seed $seed --ptype p
# done

for model in "${models[@]}";
do
    for seed in "${seeds[@]}";
    do
        printf "Config :%s\n\n" "$seed"
        python $SCRIPT_DIR/fit_model.py --source SBG_23 --detector TA2015 --model $model --dtype sim --seed $seed --ptype p --sim_model joint_gmf
    done
done

# for seed in "${seeds[@]}";
# do
#     printf "Config :%s\n\n" "$seed"
#     python $SCRIPT_DIR/fit_model.py --source 2FHL_250Mpc --detector TA2015 --model joint_gmf --dtype sim --seed $seed --ptype p --sim_model joint_gmf --tight_B
# done

# for seed in "${seeds[@]}"; 
# do
#     printf  "Config :%s\n\n" "$seed"
#     python $SCRIPT_DIR/simulate_data.py --source SBG_23 --detector TA2015 --model joint_gmf --seed $seed --ptype He
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
# for ptype in "${ptypes[@]}"; 
# do
#     for src in "${sources[@]}";
#     do
#         for seed in "${seeds[@]}";
#         do
#             printf "Config :%s\n\n" "$src joint_gmf data $ptype $seed"
#             python $SCRIPT_DIR/fit_model.py --source $src --detector TA2015 --model joint_gmf --dtype real --seed $seed --ptype $ptype
#         done
#     done
# done

# # gmf fits to data, only using joint+gmf model
# # iterate over all ptypes
# for model in "${sim_models[@]}"; 
# do
#     for ptype in "${ptypes[@]}"; 
#     do
#         for src in "${sources[@]}";
#         do
#             for seed in "${seeds[@]}";
#             do
#                 printf "Config :%s\n\n" "$src joint_gmf data $ptype $seed"
#                 python $SCRIPT_DIR/fit_model.py --source $src --detector TA2015 --model $model --dtype real --seed $seed --ptype $ptype --tight_B
#             done
#         done
#     done
# done
