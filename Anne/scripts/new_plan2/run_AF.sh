#!/usr/bin/bash

#SBATCH --job-name=AF
#SBATCH --output=AF.out
#SBATCH --error=AF.err
#SBATCH --time=160:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage


SCRIPT_PATH='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/new_plan2/'

python3 ${SCRIPT_PATH}AF_MAF_db.py