#!/usr/bin/bash
# Number or specific tissue of a sample
NUMBER=4104
# The entire file number
FILE_NUM=SS600${NUMBER}
# The path where the file is located
PATH_DICT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/S01/${NUMBER}/${METHOD}/
# Which method was used
# Name of the folder
METHOD=bwa_mem #bowtie,   bwa_aln, bwa_mem
# Piece with which the file name begins
METH_FILE=mem #bowtie2, aln, mem
# Load packages
ml SAMtools/1.9-foss-2018b

# Function that changed the sample names
change_sample_name() {
    # Change sample name
    samtools view -H ${2}${1}.bam  | sed "s/SM:[^\t]*/SM:${1}/g" | samtools reheader - ${2}${1}.bam > ${2}SN_${FILE_NUM}.bam
    # Index new file with sample names
    samtools index ${2}SN_${FILE_NUM}.bam
}
# Call function
change_sample_name ${METH_FILE}_${FILE_NUM}.DR ${PATH_DICT} #1=file
