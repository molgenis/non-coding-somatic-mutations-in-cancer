#!/bin/bash

#SBATCH --job-name=vcf_compare
#SBATCH --output=vcf_compare.out
#SBATCH --error=vcf_compare.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

# Load package
#ml BCFtools/1.11-GCCcore-7.3.0
#ml GATK/4.1.4.1-Java-8-LTS

# Loop over rows in file
while IFS= read -r line; do
    echo "tester: $line"
    # Last string in line
    echo ${line##* }
    # Compare two VCF files
    bcftools isec ${line}
    for filename in ${line##* }*.vcf; do
        bcftools view ${filename} -Oz -o ${filename}.gz
        bcftools index ${filename}.gz
    done
    echo 'END'
done < "${GENERAL_PATH}${CHROM}_compare_hc_tumor_${METHOD}_both.txt"