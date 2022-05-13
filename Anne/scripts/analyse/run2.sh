#!/usr/bin/bash

#SBATCH --job-name=region
#SBATCH --output=region.out
#SBATCH --error=region.err
#SBATCH --time=50:59:59
#SBATCH --cpus-per-task=1
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage

SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/analyse/

python3 ${SCRIPT_PATH}region.py