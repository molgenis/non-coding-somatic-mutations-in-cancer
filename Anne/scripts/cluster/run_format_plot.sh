#!/usr/bin/bash

#SBATCH --job-name=run_bowtie
#SBATCH --output=run_bowtie.out
#SBATCH --error=run_bowtie.err
#SBATCH --time=79:59:59
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
  echo "tester: $line"
  python3 ${PATH_PLOT}make_df.py ${line} ${PATH_PLOT}
done < "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/files_bowtie.txt"
#s_aln