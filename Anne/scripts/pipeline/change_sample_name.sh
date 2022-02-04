#!/usr/bin/bash

#SBATCH --job-name=change_sname
#SBATCH --output=change_sname.out
#SBATCH --error=change_sname.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

echo 'change sample name'

# Loop over TISSUE_ARR
for i in "${TISSUE_ARR[@]}"
do
    # Number or specific tissue of a sample
    NUMBER=${i}
    # The entire file number
    FILE_NUM=SS600${NUMBER}
    # The path where the file is located
    PATH_DICT=${GENERAL_PATH}${SAMPLE}/${NUMBER}/${CHROM}/${METHOD}/

    echo "BEGIN"
    echo $NUMBER
    # Function that changed the sample names
    change_sample_name() {
        samtools view -H ${2}${1}.bam  | sed "s/SM:[^\t]*/SM:${1}/g" | samtools reheader - ${2}${1}.bam > ${2}SN_${FILE_NUM}.bam
        samtools index ${2}SN_${FILE_NUM}.bam
    }

    #source /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/make_vcf.sh
    change_sample_name ${METH_FILE}_${FILE_NUM}.DR ${PATH_DICT} #1=file

    echo "EIND change sample name - ${NUMBER}"
done

