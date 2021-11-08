#!/usr/bin/bash

dir=5042_5044
TUMOR=5044_1
NORMAL=5042

PATH_DICT_one=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${TUMOR}/bwa_aln/
PATH_DICT_two=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NORMAL}/bwa_aln/
OUTPUT_DICT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NORMAL}_${TUMOR}/bwa_aln/

mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${dir}
mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${dir}/bowtie
mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${dir}/bwa_aln
mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${dir}/bwa_mem

source ./make_vcf.sh
printf "####change_sample_name\n"
#change_sample_name SS600${TUMOR}_${CHR}_mapped_1_2.bwa.RG ${PATH_DICT_one} #1=file
#change_sample_name SS600${TUMOR}_${CHR}_mapped_1_2.bwa.DR ${PATH_DICT_one} #1=file
change_sample_name bwa_aln_${NORMAL}_${CHR}_mapped_1_2.RG ${PATH_DICT_two} #1=file
change_sample_name bwa_aln_${NORMAL}_${CHR}_mapped_1_2.DR ${PATH_DICT_two} #1=file


printf "\n####Mutect2_one_sample\n"
#Mutect2_one_sample SS600${TUMOR}_${CHR}_mapped_1_2.bwa.RG tumor ${PATH_DICT_one} #1=file 2=type
#Mutect2_one_sample SS600${TUMOR}_${CHR}_mapped_1_2.bwa.DR tumor ${PATH_DICT_one} #1=file 2=type
Mutect2_one_sample bwa_aln_${NORMAL}_${CHR}_mapped_1_2.RG normal ${PATH_DICT_two} #1=file 2=type
#Mutect2_one_sample bwa_aln_${NORMAL}_${CHR}_mapped_1_2.DR normal ${PATH_DICT_two} #1=file 2=type

printf "\n####Mutect2_two_sample\n"
#Mutect2_two_samples bwa_aln_${TUMOR}_${CHR}_mapped_1_2.RG bwa_aln_${NORMAL}_${CHR}_mapped_1_2.RG ${OUTPUT_DICT} #1=tumor 2=normal
#Mutect2_two_samples SS600${TUMOR}_${CHR}_mapped_1_2.bwa.DR bwa_aln_${NORMAL}_${CHR}_mapped_1_2.DR  ${OUTPUT_DICT} #1=tumor 2=normal