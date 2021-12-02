#!/usr/bin/bash

#SBATCH --job-name=vcf_aln
#SBATCH --output=vcf_aln.out
#SBATCH --error=vcf_aln.err
#SBATCH --time=167:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


# Array of sample numbers
array2=( S1  S2  S3  S4  S5  S6 )
TYPE_ALN="mutect_bwa_aln"

CHROM=chr22

for i in "${array2[@]}"
do 
    PATH_FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/"${i}"/${CHROM}/${TYPE_ALN}/
    for filename in ${PATH_FILE}*.txt
    do
        # Name of input file
        echo ${filename}
        # Path where the file is located
        INPUT_FILE=${filename}
        # Each line in the file is an argument passed to another file.
        # Here the lines are extracted and stored as argument/variable.
        mkdir_path=$(sed '1q;d' ${INPUT_FILE})
        arg_mutect2=$(sed '2q;d' ${INPUT_FILE})
        file_output=$(sed '3q;d' ${INPUT_FILE})
        # The file from which a function is called.
        source /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/make_vcf_auto.sh
        # The function called from the file.
        # The three arguments are given.
        # The argument in quotes "" consists of multiple words (split by space), 
        # so it must be passed in quotes. Otherwise, only the first 'word' 
        # will be grabbed before the first space character.
        mutect2_vcf ${mkdir_path} "${arg_mutect2}" ${file_output}
        echo "EIND"
    done    
done

