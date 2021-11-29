#!/usr/bin/bash

#SBATCH --job-name=job_mem
#SBATCH --output=job_mem.out
#SBATCH --error=job_mem.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


# Load packages
ml SAMtools/1.9-foss-2018b
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS
ml FastQC/0.11.8-Java-11-LTS
ml cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare

# Array of tissue numbers
array=( 5042  5044 )
# Array of sample numbers
array2=( S6  S6 )

for i in "${!array[@]}"
do
    # Number or specific tissue of a sample
    NUMBER="${array[i]}"
    # The entire file number
    FILE_NUM=SS600${NUMBER}    
    # The path where the file is located
    PATH_DIR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/"${array2[i]}"/
    # The path to the file of the reference genome that will be used.
    GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa

    echo "BEGIN"
    # Creates a specific folder
    mkdir -p ${PATH_DIR}${NUMBER}/bwa_mem
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

    #BWA = Burrows-Wheeler Aligner
    #https://ucdavis-bioinformatics-training.github.io/2017-August-Variant-Analysis-Workshop/wednesday/alignment.html
    # Paired-end alignment BWA-mem
    bwa mem ${GENOOM} ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.fastq ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.fastq > ${PATH_DIR}${NUMBER}/bwa_mem/mem_${FILE_NUM}.sam
    align_last_steps ${PATH_DIR}${NUMBER}/bwa_mem/mem_${FILE_NUM}

    echo "EIND"
done

