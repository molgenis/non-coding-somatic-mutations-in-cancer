#!/bin/bash

# Load package
ml BCFtools/1.11-GCCcore-7.3.0
ml GATK/4.1.4.1-Java-8-LTS

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
    echo 'KLAAR'
done < "${GENERAL_PATH}${CHROM}_compare_hc_tumor_${METHOD}_both.txt"