#!/usr/bin/bash

#SBATCH --job-name=TFBS
#SBATCH --output=TFBS.out
#SBATCH --error=TFBS.err
#SBATCH --time=160:59:59
#SBATCH --cpus-per-task=1
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage

source /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/read_yaml.sh
YAML_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/config.yaml

SCRIPT_PATH=$(yaml $YAML_PATH "['new_plan2']")

python3 ${SCRIPT_PATH}extra_layers_TFBS.py