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

NUMBER=4104
FILE_NUM=SS600${NUMBER}
GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa
PATH_DIR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/S01/

ml SAMtools/1.9-foss-2018b
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS
ml FastQC/0.11.8-Java-11-LTS
ml cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare

echo "BEGIN"
#echo "${NUMBER}"

mkdir -p ${PATH_DIR}${NUMBER}/bwa_aln


cd ${PATH_DIR}${NUMBER}
#echo "BEGIN sort"
# BAM files must be resorted so that they are ordered by read ID instead of location in the reference
#samtools view -h ${PATH_DIR}${FILE_NUM}.sorted.bam chr22 > ${PATH_DIR}${FILE_NUM}_chr22.sam
#samtools view -bS ${PATH_DIR}${FILE_NUM}_chr22.sam > ${PATH_DIR}${FILE_NUM}_chr22.bam
#samtools sort -n ${PATH_DIR}${FILE_NUM}_chr22.bam -o ${PATH_DIR}${NUMBER}/${FILE_NUM}_byName.sorted.bam
# extract the FASTQ reads into two paired read files
#echo "KLAAR sort"
#echo "BEGIN 1 en 2"
#cd ${PATH_DIR}${NUMBER}
#samtools fastq -@ 8 ${PATH_DIR}${NUMBER}/${FILE_NUM}_byName.sorted.bam -1 ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.fastq -2 ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.fastq -0 /dev/null -s /dev/null -n

#fastqc -f fastq -o ${PATH_DIR}${NUMBER}/QC ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.fastq ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.fastq
#fastqc -f bam -o ${PATH_DIR}${NUMBER}/QC ${PATH_DIR}${NUMBER}/${FILE_NUM}_byName.sorted.bam

#echo "KLAAR 1 en 2"

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



#echo "BWA aln"
bwa aln ${GENOOM} ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.fastq > ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.sai && bwa aln ${GENOOM} ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.fastq > ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.sai && bwa sampe ${GENOOM} ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.sai ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.sai ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R1.fastq ${PATH_DIR}${NUMBER}/${FILE_NUM}_name_R2.fastq > ${PATH_DIR}${NUMBER}/bwa_aln/aln_${FILE_NUM}.sam
align_last_steps ${PATH_DIR}${NUMBER}/bwa_aln/aln_${FILE_NUM}


echo "EIND"