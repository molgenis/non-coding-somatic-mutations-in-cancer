#!/usr/bin/bash

dir=5042_5044
TUMOR=5044_1
NORMAL=5042

PATH_DICT_one=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${TUMOR}/bowtie/
PATH_DICT_two=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NORMAL}/bowtie/
OUTPUT_DICT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NORMAL}_${TUMOR}/bowtie/

mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${dir}
mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${dir}/bowtie
mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${dir}/bwa_aln
mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${dir}/bwa_mem

source ./make_vcf.sh
printf "####change_sample_name\n"
change_sample_name bowtie2_aln_pe.RG ${PATH_DICT_one} #1=file
change_sample_name bowtie2_aln_pe.DR ${PATH_DICT_one} #1=file
change_sample_name bowtie2_aln_pe.RG ${PATH_DICT_two} #1=file
change_sample_name bowtie2_aln_pe.DR ${PATH_DICT_two} #1=file

printf "\n####Mutect2_one_sample\n"
#Mutect2_one_sample bowtie2_aln_pe.RG tumor ${PATH_DICT_one} #1=file 2=type
#Mutect2_one_sample bowtie2_aln_pe.DR tumor ${PATH_DICT_one} #1=file 2=type
Mutect2_one_sample bowtie2_aln_pe.RG normal ${PATH_DICT_two} #1=file 2=type
Mutect2_one_sample bowtie2_aln_pe.DR normal ${PATH_DICT_two} #1=file 2=type

printf "\n####Mutect2_two_sample\n"
Mutect2_two_samples bowtie2_aln_pe.RG bowtie2_aln_pe.RG ${OUTPUT_DICT} #1=tumor 2=normal
Mutect2_two_samples bowtie2_aln_pe.DR bowtie2_aln_pe.DR  ${OUTPUT_DICT} #1=tumor 2=normal