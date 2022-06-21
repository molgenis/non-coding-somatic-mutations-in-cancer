#!/usr/bin/bash

#SBATCH --job-name=DNase2
#SBATCH --output=DNase2.out
#SBATCH --error=DNase2.err
#SBATCH --time=160:59:59
#SBATCH --cpus-per-task=16
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage

SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/add_layers/

python3 ${SCRIPT_PATH}extra_layers_DNase_multi.py