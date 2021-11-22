#!/usr/bin/bash

# Load packages
ml BCFtools/1.11-GCCcore-7.3.0

PATH_GENERAL=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/merge_vcf/
OUTPUT_PATH=${PATH_GENERAL}annotated/
#bcftools annotate -a /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/dbSNP/merge_All_20180423.vcf.gz   -o /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/00112test_5044_OUTPUT.vcf  /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/00test_5044.vcf.gz  

# Loop over all .csv files in this folder
for filename in ${PATH_GENERAL}/*.vcf; do
    echo ${filename}
    echo "$(basename -- $filename)"
    # BASENAME=$( echo "$(basename -- $filename)")
    # bcftools view ${filename} -Oz -o ${filename}.gz
    # bcftools index ${filename}.gz
    # #-c CHROM,FROM,TO,ID 
    # bcftools annotate -c CHROM,FROM,TO,ID -a /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/dbSNP/merge_All_20180423.vcf.gz   -o ${OUTPUT_PATH}m_${BASENAME}  ${filename}.gz

    bcftools view --header-only ${OUTPUT_PATH}m_${BASENAME} > ${OUTPUT_PATH}dbSNP_filter/header_${BASENAME}
    python3 filter_dbSNP.py ${OUTPUT_PATH}m_${BASENAME} ${OUTPUT_PATH}dbSNP_filter/
    cat ${OUTPUT_PATH}dbSNP_filter/header_${BASENAME} ${OUTPUT_PATH}dbSNP_filter/noHeader_${BASENAME} >> ${OUTPUT_PATH}dbSNP_filter/merge_${BASENAME}
done

echo 'END'