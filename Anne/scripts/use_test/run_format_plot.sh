#!/usr/bin/bash

#SBATCH --job-name=run_aln
#SBATCH --output=run_aln.out
#SBATCH --error=run_aln.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

# Load packages
ml Anaconda3/5.3.0
source activate stage

PATH_GLOBAL=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/
PATH_PLOT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/FORMAT/

# Loop over lines in file
while IFS= read -r line; do
  # Call python file
  # The python file creates dataframe from the manually created tumor files 
  # (by comparing hc against tumor, and grabbing the unique tumor SNPs)
  python3 ${PATH_PLOT}make_df.py ${line} ${PATH_PLOT}
done < "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/files_aln.txt"