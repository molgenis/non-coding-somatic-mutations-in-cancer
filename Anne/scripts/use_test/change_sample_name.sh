#!/usr/bin/bash

NUMBER=4104
FILE_NUM=SS600${NUMBER}
METHOD=bwa_mem #bowtie,   bwa_aln, bwa_mem
METH_FILE=mem #bowtie2, aln, mem
PATH_DICT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/S01/${NUMBER}/${METHOD}/

ml SAMtools/1.9-foss-2018b

# Function that changed the sample names
change_sample_name() {
    samtools view -H ${2}${1}.bam  | sed "s/SM:[^\t]*/SM:${1}/g" | samtools reheader - ${2}${1}.bam > ${2}SN_${FILE_NUM}.bam
    samtools index ${2}SN_${FILE_NUM}.bam
}

#source /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/make_vcf.sh
#printf "####change_sample_name\n"
change_sample_name ${METH_FILE}_${FILE_NUM}.DR ${PATH_DICT} #1=file

#bowtie2_SS6004094.DR.bam