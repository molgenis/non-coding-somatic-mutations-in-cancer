#!/bin/bash

# Then each VCF file of a tumor is compared
# with that of the corresponding germline variant using BCFtools. 
# This results in three files: 1) with unique SNPs from the tumor, 
#                              2) with unique SNPs from the germline and 
#                              3) with SNPs shared by both (overlap). 
# This was done with the goal of seeing if we got the same results
# using a tumor filtered by the germline variation, as we did by running it with one tumor and one
# germline variant together. 

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

# Loop over rows in file
while IFS= read -r line; do
    echo "tester: $line"
    # Last string in line
    echo ${line##* }
    # Compare two VCF files
    bcftools isec ${line}
    # Loop over lines that ends with .vcf
    for filename in ${line##* }*.vcf; do
        echo ${filename}
        echo 'COMPARE'
        # Make .gz file of file
        bcftools view ${filename} -Oz -o ${filename}.gz
        # Index file
        bcftools index ${filename}.gz
    done
    echo 'END'
done < "${GENERAL_PATH}${CHROM}_compare_hc_tumor_${METHOD}_both.txt"