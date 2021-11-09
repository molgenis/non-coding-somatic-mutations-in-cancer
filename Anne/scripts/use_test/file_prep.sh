#!/usr/bin/bash

# Load packages
# The path to the file of the reference genome that will be used.
GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa
# Load packages
ml SAMtools/1.9-foss-2018b
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS
ml FastQC/0.11.8-Java-11-LTS
ml cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare

# Array of files it should download
array=( 4109  4113  4114  4118  4119  4123  4124  4129  4128  4133  4134  4139  4138  5041  5043  5042  5044 )
# Array with file output names
array2=( S2  S2  S2  S2  S2  S3  S3  S3  S4  S4  S4  S5  S5  S5  S6  S6  S6 )

METHOD=bwa_aln #bowtie,   bwa_aln, bwa_mem
METH_FILE=aln #bowtie2, aln, mem

for i in "${!array[@]}"
do
    # Number or specific tissue of a sample
    NUMBER="${array[i]}"
    # The entire file number
    FILE_NUM=SS600${NUMBER}
    # The path where the file is located
    PATH_DIR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/"${array2[i]}"/

    echo "BEGIN"
    # Creates a specific folder
    mkdir -p ${PATH_DIR}${NUMBER}/QC
    cd ${PATH_DIR}${NUMBER}
    # Filters a certain chromosome from the bam file
    samtools view -h ${PATH_DIR}${FILE_NUM}.sorted.bam chr22 > ${PATH_DIR}${FILE_NUM}_chr22.sam
    samtools view -bS ${PATH_DIR}${FILE_NUM}_chr22.sam > ${PATH_DIR}${FILE_NUM}_chr22.bam
    # BAM files must be resorted so that they are ordered by read ID instead of location in the reference
    samtools sort -n ${PATH_DIR}${FILE_NUM}_chr22.bam -o ${PATH_DIR}${NUMBER}/${FILE_NUM}_byName.sorted.bam
    # Extract the FASTQ reads into two paired read files
    cd ${PATH_DIR}${NUMBER}
    samtools fastq -@ 8 ${PATH_DIR}${NUMBER}/${FILE_NUM}_byName.sorted.bam -1 ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.fastq -2 ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.fastq -0 /dev/null -s /dev/null -n
    # Runs FastQC to see the quality of the reads of the file.
    # On fastq files
    fastqc -f fastq -o ${PATH_DIR}${NUMBER}/QC ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.fastq ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.fastq
    # On bam file
    fastqc -f bam -o ${PATH_DIR}${NUMBER}/QC ${PATH_DIR}${NUMBER}/${FILE_NUM}_byName.sorted.bam

    echo "EIND"
done
