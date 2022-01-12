#!/usr/bin/bash

# Load packages
ml BCFtools/1.11-GCCcore-7.3.0


#seq FIRST STEP LAST
chrom_num=($(seq 1 1 22))
chrom_num+=("X" "Y")

PATH_DBSNP=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/dbSNP/

for i in "${chrom_num[@]}"
do
    #for ((i=1;i<=22;i++)); do
    echo $i
    bcftools view ${PATH_DBSNP}merge_All_20180423.vcf.gz --regions chr${i} > ${PATH_DBSNP}per_chr/chr${i}_merge_All_20180423.vcf
    #gzip -c ${PATH_DBSNP}per_chr/chr${i}_merge_All_20180423.vcf > ${PATH_DBSNP}per_chr/chr${i}_merge_All_20180423.vcf.gz
    bgzip ${PATH_DBSNP}per_chr/chr${i}_merge_All_20180423.vcf #bcftools view file.vcf -Oz -o file.vcf.gz
    tabix ${PATH_DBSNP}per_chr/chr${i}_merge_All_20180423.vcf.gz #bcftools index file.vcf.gz

done



echo 'END'