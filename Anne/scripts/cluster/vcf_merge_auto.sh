#!/bin/bash

#SBATCH --job-name=vcf_merge
#SBATCH --output=vcf_merge.out
#SBATCH --error=vcf_merge.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

# Path of the file that will be used
FILE_PATH=${GENERAL_PATH}${COMP_TYPE}_comparison_${METHOD}_both.txt
# 
MERGE_COM=$(sed '1q;d' ${FILE_PATH})
echo ${MERGE_COM}
mkdir -p ${GENERAL_PATH}merge_vcf/${CHROM}/
#--force-samples: if the merged files contain duplicate samples names, proceed anyway. 
# Duplicate sample names will be resolved by prepending the index of the file as it 
# appeared on the command line to the conflicting sample name (see 2:S3 in the above example).
bcftools merge ${MERGE_COM} -o ${GENERAL_PATH}merge_vcf/${CHROM}/${COMP_TYPE}_${METHOD}.vcf --force-samples