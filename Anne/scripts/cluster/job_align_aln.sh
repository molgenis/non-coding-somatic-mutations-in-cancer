#!/usr/bin/bash

#SBATCH --job-name=job_aln
#SBATCH --output=job_aln.out
#SBATCH --error=job_aln.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L




for i in "${!array[@]}"
do
    # Number or specific tissue of a sample
    NUMBER="${array[i]}"
    # The entire file number
    FILE_NUM=SS600${NUMBER}  
    # The path where the file is located
    PATH_DIR=${GENERAL_PATH}"${array2[i]}"/
    # The path to the file of the reference genome that will be used.
    #GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHROM}.fa #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa  

    echo "BEGIN"
    echo $NUMBER
    # Creates a specific folder
    mkdir -p ${PATH_DIR}${NUMBER}/${CHROM}/bwa_aln/
    cd ${PATH_DIR}${NUMBER}

    # Last steps of alignment
    align_last_steps() {
        samtools view -Sb ${1}.sam > ${1}.bam
        print('TEST1')
        samtools sort ${1}.bam -o ${1}_sort.bam
        print('TEST2')
        samtools index ${1}_sort.bam
        print('TEST3')
        # Adding Read Group tags and indexing bam files
        # ml picard/2.20.5-Java-11-LTS
        print('TEST4')
        java -jar /apps/software/picard/2.20.5-Java-11-LTS/picard.jar  AddOrReplaceReadGroups INPUT= ${1}_sort.bam OUTPUT= ${1}.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT
        samtools index ${1}.RG.bam
        print('TEST5')
        # Marking and removing duplicates
        java -jar /apps/software/picard/2.20.5-Java-11-LTS/picard.jar  MarkDuplicates I= ${1}.RG.bam O= ${1}.DR.bam M=${1}_output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> ${1}_logFile.log
        samtools index ${1}.DR.bam
    }

    # BWA = Burrows-Wheeler Aligner
    # Paired-end alignment BWA-aln
    bwa aln ${GENOOM} ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.fastq > ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.sai && bwa aln ${GENOOM} ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.fastq > ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.sai && bwa sampe ${GENOOM} ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.sai ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.sai ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.fastq ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.fastq > ${PATH_DIR}${NUMBER}/${CHROM}/bwa_aln/aln_${FILE_NUM}.sam
    print('begin alignen')
    align_last_steps ${PATH_DIR}${NUMBER}/${CHROM}/bwa_aln/aln_${FILE_NUM}

    echo "EIND"
done

