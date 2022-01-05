#!/usr/bin/bash

#SBATCH --job-name=file_prep
#SBATCH --output=file_prep.out
#SBATCH --error=file_prep.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

echo 'file prep'

for i in "${TISSUE_ARR[@]}"
do
    # Number or specific tissue of a sample
    NUMBER=${i}
    # The entire file number
    FILE_NUM=SS600${NUMBER}
    # The path where the file is located
    PATH_DIR=${GENERAL_PATH}${SAMPLE}/

    echo "BEGIN"
    echo ${NUMBER}
    # Creates a specific folder
    mkdir -p ${PATH_DIR}${NUMBER}/${CHROM}/QC

    cd ${PATH_DIR}${NUMBER}
    # Filters a certain chromosome from the bam file
    samtools view -h ${PATH_DIR}${FILE_NUM}.sorted.bam ${CHROM} > ${PATH_DIR}${FILE_NUM}_${CHROM}.sam
    samtools view -bS ${PATH_DIR}${FILE_NUM}_${CHROM}.sam > ${PATH_DIR}${FILE_NUM}_${CHROM}.bam
    # BAM files must be resorted so that they are ordered by read ID instead of location in the reference
    samtools sort -n ${PATH_DIR}${FILE_NUM}_${CHROM}.bam -o ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_byName.sorted.bam
    # Extract the FASTQ reads into two paired read files
    cd ${PATH_DIR}${NUMBER}
    samtools fastq -@ 8 ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_byName.sorted.bam -1 ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.fastq -2 ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.fastq -0 /dev/null -s /dev/null -n
    # Runs FastQC to see the quality of the reads of the file.
    # On fastq files
    fastqc -f fastq -o ${PATH_DIR}${NUMBER}/${CHROM}/QC ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.fastq ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.fastq
    # On bam file
    fastqc -f bam -o ${PATH_DIR}${NUMBER}/${CHROM}/QC ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_byName.sorted.bam

    echo "EIND file prep - ${NUMBER}"
done

