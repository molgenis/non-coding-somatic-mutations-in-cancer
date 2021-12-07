#!/usr/bin/bash

#SBATCH --job-name=job_bowtie
#SBATCH --output=job_bowtie.out
#SBATCH --error=job_bowtie.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L



#METHOD=bwa_aln #bowtie,   bwa_aln, bwa_mem
#METH_FILE=aln #bowtie2, aln, mem

for i in "${!array[@]}"
do
    # Number or specific tissue of a sample
    NUMBER="${array[i]}"
    echo ${NUMBER}
    FILE_NUM=SS600${NUMBER}
    PATH_DICT=${GENERAL_PATH}"${array2[i]}"/${NUMBER}/${CHROM}/${METHOD}/

    # Function that changed the sample names
    change_sample_name() {
        samtools view -H ${2}${1}.bam  | sed "s/SM:[^\t]*/SM:${1}/g" | samtools reheader - ${2}${1}.bam > ${2}SN_${FILE_NUM}.bam
        samtools index ${2}SN_${FILE_NUM}.bam
    }

    #source /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/make_vcf.sh
    #printf "####change_sample_name\n"
    change_sample_name ${METH_FILE}_${FILE_NUM}.DR ${PATH_DICT} #1=file

    #bowtie2_SS6004094.DR.bam

    echo "EIND"
done

