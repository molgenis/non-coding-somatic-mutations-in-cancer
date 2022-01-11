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

PATH_PLOT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/FORMAT/

# Loop over all .csv files in this folder
for filename in ${PATH_PLOT}/*.csv; do
    # Get file path without last extension (.csv)
    NEW_FOLDER="${filename%.*}"
    # Create new folder, to save all images in
    mkdir -p ${NEW_FOLDER}
    # Call python file
    # Python file creates histograms of the entered data per column
    python3 ${PATH_PLOT}make_plots.py ${filename} ${NEW_FOLDER}/
done