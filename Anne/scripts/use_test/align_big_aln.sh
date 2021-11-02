#!/usr/bin/bash

#SBATCH --job-name=align_big
#SBATCH --output=align_big.out
#SBATCH --error=align_big.err
#SBATCH --time=167:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

# Number or specific tissue of a sample
NUMBER=4104
# The entire file number
FILE_NUM=SS600${NUMBER}
# The path where the file is located
PATH_DIR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/S01/
# The path to the file of the reference genome that will be used.
GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa
# Load packages
ml SAMtools/1.9-foss-2018b
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS
ml FastQC/0.11.8-Java-11-LTS
ml cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare

echo "BEGIN"
# Creates a specific folder
mkdir -p ${PATH_DIR}${NUMBER}/bwa_aln
cd ${PATH_DIR}${NUMBER}

# Last steps of alignment
align_last_steps() {
    samtools view -Sb ${1}.sam > ${1}.bam
    samtools sort ${1}.bam -o ${1}_sort.bam
    samtools index ${1}_sort.bam
    # Adding Read Group tags and indexing bam files
    java -jar ${EBROOTPICARD}/picard.jar  AddOrReplaceReadGroups INPUT= ${1}_sort.bam OUTPUT= ${1}.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT
    samtools index ${1}.RG.bam
    # Marking and removing duplicates
    java -jar ${EBROOTPICARD}/picard.jar  MarkDuplicates I= ${1}.RG.bam O= ${1}.DR.bam M=${1}_output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> ${1}_logFile.log
    samtools index ${1}.DR.bam
}

# BWA = Burrows-Wheeler Aligner
# Paired-end alignment BWA-aln
bwa aln ${GENOOM} ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.fastq > ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.sai && bwa aln ${GENOOM} ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.fastq > ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.sai && bwa sampe ${GENOOM} ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.sai ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.sai ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.fastq ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.fastq > ${PATH_DIR}${NUMBER}/bwa_aln/aln_${FILE_NUM}.sam
align_last_steps ${PATH_DIR}${NUMBER}/bwa_aln/aln_${FILE_NUM}

echo "EIND"