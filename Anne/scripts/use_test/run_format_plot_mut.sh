#!/usr/bin/bash

#SBATCH --job-name=run_bowtie_mut
#SBATCH --output=run_bowtie_mut.out
#SBATCH --error=run_bowtie_mut.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage

PATH_GLOBAL=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/
PATH_PLOT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/FORMAT/


while IFS= read -r line; do
  echo ${line}
  echo ""
  NEW_FILE=${line%.*}
  echo ${NEW_FILE}
  echo ""
  zcat "${NEW_FILE}".gz > "${NEW_FILE}"
  python3 ${PATH_PLOT}make_df_mut.py ${NEW_FILE} ${PATH_PLOT}mut_
done < "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/files_mut_bowtie.txt"
#s_mut_bowtie