#!/bin/bash

# directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Initialization of arrays
sources=('SBG_23' '2FHL_250Mpc' 'swift_BAT_213')
detectors=('TA2015' 'auger2014')
models=('arrival_direction' 'joint' 'joint_gmf')
sim_models=('joint' 'joint_gmf')
dtypes=('sim' 'real')
seeds=('19990308' '16852056' '65492186' '09953497' '9999999')
ptypes=('p' 'He' 'N' 'O' 'Si' 'Fe')

# # precompute tables
# printf "\nPrecomputing tables...\n"
# for src in "${sources[@]}"; 
# do
#     printf "Current source: %s\n\n" "$src"
#     for detector in "${detectors[@]}"; 
#     do
#         printf "Current detector: %s\n\n" "$detector"
#         python $SCRIPT_DIR/tables/precompute_tables.py --source $src --detector $detector
#     done
# done

# # simulate data for all sources / detectors / models / seeds
# printf "\n Simulating UHECRs...\n"
# for sim_model in "${sim_models[@]}"; 
# do
#     printf "Current Model: %s\n\n" "$sim_model"
#     for src in "${sources[@]}"; 
#     do
#         printf "Current source: %s\n\n" "$src"
#         for detector in "${detectors[@]}"; 
#         do
#             printf "Current detector: %s\n\n" "$detector"
#             for seed in "${seeds[@]}"; 
#             do
#                 printf "Current seed: %s\n" "$seed"
#                 python $SCRIPT_DIR/simulate_data.py --source $src --detector $detector --model $sim_model --seed $seed
#             done
#             printf "\n"
#         done
#     done
# done

for seed in "${seeds[@]}"; 
do
    printf "Current seed: %s\n" "$seed"
    python $SCRIPT_DIR/fit_model.py --source "swift_BAT_213" --detector "auger2014" --model "joint" --dtype "real" --seed $seed
done


# # fit for all sources / detectors / models / seeds / dtype
# printf "\nFitting model to data...\n"
# for model in "${models[@]}"; 
# do
#     printf "Current Model: %s\n\n" "$model"
#     for src in "${sources[@]}"; 
#     do
#         printf "Current source: %s\n\n" "$src"
#         for detector in "${detectors[@]}"; 
#         do
#             printf "Current detector: %s\n\n" "$detector"
#             for dtype in "${dtypes[@]}"; 
#             do
#                 printf "Current data type used: %s\n\n" "$dtype"
#                 for seed in "${seeds[@]}"; 
#                 do
#                     printf "Current seed: %s\n" "$seed"
#                     # Run compositional arguments only when joint+gmf and real data (gmf_sims dont contain compositional info)
#                     if [[ $model == 'joint_gmf' ]] && [[ $dtype == 'real' ]];
#                     then
#                         for ptype in "${ptypes[@]}";
#                         do
#                             printf "Current particle type: %s\n" "$ptype"
#                             python $SCRIPT_DIR/fit_model.py --source $src --detector $detector --model $model --dtype $dtype --seed $seed --ptype $ptype
#                         done
#                     else
#                         printf "Current particle type: %s\n" "p"
#                         python $SCRIPT_DIR/fit_model.py --source $src --detector $detector --model $model --dtype $dtype --seed $seed
#                     fi
#                 done
#                 printf "\n"
#             done
#         done
#     done
# done


# do all gmf-related sim + fitting

# gmf sims
# for src in "${sources[@]}"; 
# do
#     printf "Current source: %s\n\n" "$src"
#     for seed in "${seeds[@]}"; 
#     do
#         printf "Current seed: %s\n" "$seed"
#         python $SCRIPT_DIR/simulate_data.py --source $src --detector TA2015 --model joint_gmf --seed $seed
#     done
#     printf "\n"
# done




# gmf fits
# for src in "${sources[@]}"; 
# do
#     for ptype in "${ptypes[@]}";
#     do
        # for dtype in "${dtypes[@]}"; 
        # do
            # for model in "${models[@]}"; 
            # do
                # for seed in "${seeds[@]}";
                # do
                    # printf "Config :%s\n\n" "$src $dtype $model $dtype $seed"
                #     python $SCRIPT_DIR/fit_model.py --source $src --detector TA2015 --model $model --sim_model joint_gmf --dtype $dtype --seed $seed --ptype $ptype
                # done
            # done
        # done
    # done
# done