#!/bin/bash

# Load package
ml BCFtools/1.11-GCCcore-7.3.0

COMP_TYPE=mutect2 # manual, mutect2
METH_TYPE=bwa_mem # bowtie, bwa_aln, bwa_mem
PATH_GL=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/
FILE_PATH=${PATH_GL}${COMP_TYPE}_comparison_${METH_TYPE}_both.txt

MERGE_COM=$(sed '1q;d' ${FILE_PATH})
echo ${MERGE_COM}
#--force-samples: if the merged files contain duplicate samples names, proceed anyway. Duplicate sample names will be resolved by prepending the index of the file as it appeared on the command line to the conflicting sample name (see 2:S3 in the above example).
bcftools merge ${MERGE_COM} -o ${PATH_GL}merge_vcf/${COMP_TYPE}_${METH_TYPE}.vcf --force-samples