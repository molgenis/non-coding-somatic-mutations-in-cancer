#!/usr/bin/bash

#SBATCH --job-name=2tables
#SBATCH --output=2tables.out
#SBATCH --error=2tables.err
#SBATCH --time=80:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=116gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage


SCRIPT_PATH='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/new_plan/'

python3 ${SCRIPT_PATH}allel_freq.py
echo 'END END END'