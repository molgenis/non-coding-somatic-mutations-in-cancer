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

#PATH_PLOT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/FORMAT/

#for filename in ${PATH_PLOT}/*.csv; do
#    NEW_FOLDER="${filename%.*}"
#    echo ${NEW_FOLDER}
#    mkdir -p ${NEW_FOLDER}
#    echo ${filename}
#    python3 ${PATH_PLOT}make_plots.py ${filename} ${NEW_FOLDER}/
#done


while IFS= read -r line; do
  OUTPUT_FILE=$(echo $line | rev | cut -d ' ' -f '1' | rev)
  mkdir -p "${OUTPUT_FILE%.*}"
  python3 ${PATH_PLOT}make_plots.py ${OUTPUT_FILE} "${OUTPUT_FILE%.*}"/
done < "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/FORMAT/part_command.txt"

echo 'END'