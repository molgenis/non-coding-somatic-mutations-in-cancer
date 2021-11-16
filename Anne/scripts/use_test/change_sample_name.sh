#!/usr/bin/bash

# Load packages
ml SAMtools/1.9-foss-2018b

# Array of tissue numbers
array=( 4109  4113  4114  4118  4119  4123  4124  4129 )
# Array of sample numbers
array2=( S2  S2  S2  S2  S2  S3  S3  S3)

METHOD=bwa_aln #bowtie,   bwa_aln, bwa_mem
METH_FILE=aln #bowtie2, aln, mem

for i in "${!array[@]}"
do
    # Number or specific tissue of a sample
    NUMBER="${array[i]}"
    echo ${NUMBER}
    # The entire file number
    FILE_NUM=SS600${NUMBER}
    # The path where the file is located
    PATH_DICT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/"${array2[i]}"/${NUMBER}/${METHOD}/
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
    echo "EIND"
done

