#!/bin/bash

FILE_PATH=${GENERAL_PATH}${COMP_TYPE}_comparison_${METHOD}_both.txt
MERGE_COM=$(sed '1q;d' ${FILE_PATH})
echo ${MERGE_COM}
#--force-samples: if the merged files contain duplicate samples names, proceed anyway. Duplicate sample names will be resolved by prepending the index of the file as it appeared on the command line to the conflicting sample name (see 2:S3 in the above example).
bcftools merge ${MERGE_COM} -o ${GENERAL_PATH}merge_vcf/${CHROM}/${COMP_TYPE}_${METHOD}.vcf --force-samples