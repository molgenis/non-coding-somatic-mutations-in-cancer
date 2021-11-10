#!/usr/bin/bash

# Load package
ml BCFtools/1.11-GCCcore-7.3.0

PATH_GL=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/

bcftools merge ${PATH_GL}S1/compair_4099_4104/bowtie/0002.vcf.gz ${PATH_GL}S1/compair_4099_4104/bowtie/0003.vcf.gz -o ${PATH_GL}merge_vcf.vcf
#https://www.biostars.org/p/392994/
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%TGT]\n' ${PATH_GL}merge_vcf.vcf | awk -v FS="\t" -v OFS="\t" '{for(i=6;i<=NF;i++) {split($i, gt, "/"); if(gt[1]==".") $i="-"; else if(gt[1]==gt[2]) $i=gt[1]; else $i="N";} print }' > ${PATH_GL}merge_vcf.txt