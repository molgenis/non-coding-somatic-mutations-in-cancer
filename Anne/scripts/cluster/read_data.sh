#!/usr/bin/bash

#SBATCH --job-name=vcf_aln
#SBATCH --output=vcf_aln.out
#SBATCH --error=vcf_aln.err
#SBATCH --time=167:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


# Array with file output names
#array2=( S2  S3  S4  S5  S6 )


PATH_FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/
SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/Anne/non-coding-somatic-mutations-in-cancer/Anne/scripts/cluster/
for filename in ${PATH_FILE}download*.tsv #.gz
do
    # Name of input file
    #echo ${filename}
    #gunzip ${filename}
    ml Anaconda3/5.3.0
    source activate stage
    python3 ${SCRIPT_PATH}read_cancer_data.py ${filename} ${PATH_FILE}
    #echo ${filename%/*}
    #echo "$(basename -- $filename)"
    echo "EIND"
done    

