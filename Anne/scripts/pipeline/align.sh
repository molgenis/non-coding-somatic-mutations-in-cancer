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

echo 'job align aln'

# Loop over TISSUE_ARR
for i in "${TISSUE_ARR[@]}"
do
    # Number or specific tissue of a sample
    NUMBER=${i}
    # The entire file number
    FILE_NUM=SS600${NUMBER}  
    # The path where the file is located
    PATH_DIR=${GENERAL_PATH}${SAMPLE}/

    echo "BEGIN"
    echo $NUMBER
    # Creates a specific folder
    mkdir -p ${PATH_DIR}${NUMBER}/${CHROM}/${METHOD}/
    cd ${PATH_DIR}${NUMBER}

    # Last steps of alignment
    align_last_steps() {
        # Load samtools (again otherwise error)
        ml SAMtools/1.9-foss-2018b
        # Make sam file into bam file 
        samtools view -Sb ${1}.sam > ${1}.bam
        rm ${1}.sam
        # Sort bam file
        samtools sort ${1}.bam -o ${1}_sort.bam
        rm ${1}.bam
        # Index bam file
        samtools index ${1}_sort.bam
        # Adding Read Group tags and indexing bam files
        # ml picard/2.20.5-Java-11-LTS
        java -jar /apps/software/picard/2.20.5-Java-11-LTS/picard.jar  AddOrReplaceReadGroups INPUT= ${1}_sort.bam OUTPUT= ${1}.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT
        rm ${1}_sort.bam
        # Index bam file
        samtools index ${1}.RG.bam
        # Marking and removing duplicates
        java -jar /apps/software/picard/2.20.5-Java-11-LTS/picard.jar  MarkDuplicates I= ${1}.RG.bam O= ${1}.DR.bam M=${1}_output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> ${1}_logFile.log
        rm ${1}.RG.bam
        rm ${1}.RG.bam.bai
        # Index bam file
        samtools index ${1}.DR.bam
    }

    if [ "${METHOD}" == "bwa_aln" ]; then
        # BWA = Burrows-Wheeler Aligner
        # Paired-end alignment BWA-aln
        bwa aln ${GENOOM} ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.fastq > ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.sai && bwa aln ${GENOOM} ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.fastq > ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.sai && bwa sampe ${GENOOM} ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.sai ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.sai ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.fastq ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.fastq > ${PATH_DIR}${NUMBER}/${CHROM}/bwa_aln/aln_${FILE_NUM}.sam
        echo 'Begin alignen'
        align_last_steps ${PATH_DIR}${NUMBER}/${CHROM}/${METHOD}/aln_${FILE_NUM}
        # Remove files
        rm ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.sai
        rm ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.sai
    elif [ "${METHOD}" == "bwa_mem" ]; then
        #BWA = Burrows-Wheeler Aligner
        #https://ucdavis-bioinformatics-training.github.io/2017-August-Variant-Analysis-Workshop/wednesday/alignment.html
        # Paired-end alignment BWA-mem
        bwa mem ${GENOOM} ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.fastq ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.fastq > ${PATH_DIR}${NUMBER}/${CHROM}/${METHOD}/mem_${FILE_NUM}.sam
        align_last_steps ${PATH_DIR}${NUMBER}/${CHROM}/${METHOD}/mem_${FILE_NUM}   
    elif [ "${METHOD}" == "bowtie" ]; then
        #BOWTIE2
        # Paired-end alignment bowtie2
        bowtie2-build ${GENOOM} ${GENERAL_PATH}GENOOM_${CHROM}
        bowtie2 -p 4 -x ${GENERAL_PATH}GENOOM_${CHROM} -1 ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.fastq -2 ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.fastq -S ${PATH_DIR}${NUMBER}/${CHROM}/${METHOD}/bowtie2_${FILE_NUM}.sam
        align_last_steps ${PATH_DIR}${NUMBER}/${CHROM}/${METHOD}/bowtie2_${FILE_NUM}
    else
        echo "ERROR"
    fi

    # Remove files
    rm ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R1.fastq
    rm ${PATH_DIR}${NUMBER}/${CHROM}/${FILE_NUM}_name_R2.fastq

    echo "EIND job align ${METHOD} - ${NUMBER}"
done

