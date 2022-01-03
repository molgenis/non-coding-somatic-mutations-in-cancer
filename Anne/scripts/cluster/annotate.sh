#!/usr/bin/bash

#SBATCH --job-name=annotate
#SBATCH --output=annotate.out
#SBATCH --error=annotate.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

# Load packages
#ml BCFtools/1.11-GCCcore-7.3.0

PATH_GENERAL=${GENERAL_PATH}merge_vcf/${CHROM}/
OUTPUT_PATH=${PATH_GENERAL}annotated/

mkdir -p ${PATH_GENERAL}annotated/
mkdir -p ${PATH_GENERAL}dbSNP_filter/
#bcftools annotate -a /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/dbSNP/merge_All_20180423.vcf.gz   -o /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/00112test_5044_OUTPUT.vcf  /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/00test_5044.vcf.gz  

# Loop over all .vcf files in this folder
for filename in ${PATH_GENERAL}/*.vcf; do
    echo ${filename}
    echo "$(basename -- $filename)"
    # Get basename of file
    BASENAME=$( echo "$(basename -- $filename)")
    bcftools view ${filename} -Oz -o ${filename}.gz
    bcftools index ${filename}.gz
    #-c CHROM,FROM,TO,ID 
    bcftools annotate -c CHROM,FROM,TO,ID -a /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/dbSNP/merge_All_20180423.vcf.gz   -o ${OUTPUT_PATH}${BASENAME}  ${filename}.gz
    # Grab the header
    bcftools view --header-only ${OUTPUT_PATH}${BASENAME} > ${PATH_GENERAL}dbSNP_filter/header_${BASENAME}
    # Filter the vcf file (with dbSNP)
    python3 ${PATH_GENERAL}filter_dbSNP.py ${OUTPUT_PATH}${BASENAME} ${PATH_GENERAL}dbSNP_filter/
    # Merge header and noHeader
    cat ${PATH_GENERAL}dbSNP_filter/header_${BASENAME} ${PATH_GENERAL}dbSNP_filter/noHeader_${BASENAME} >> ${PATH_GENERAL}dbSNP_filter/merge_${BASENAME}
done

echo 'END'