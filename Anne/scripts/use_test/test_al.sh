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

# The path where the file is located
PATH_DIR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/"${array2[i]}"/
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
array=( 4099  4094  4104  4113  4114  4118  4109  4119  4124  4123  4129  4133  4128  4134  4139
4138  5041  5043  5042  5044 )
# Array with file output names
array2=( S1  S1  S1  S2  S2  S2  S2  S2  S3  S3  S3  S4  S4  S4  S5
S5  S5  S6  S6  S6 )
TYPE_ALN="mutect_bowtie"

for i in "${!array[@]}"
do 
    PATH_FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/"${array2[i]}"/${TYPE_ALN}/
    for filename in ${PATH_FILE}*.txt
    do
        # Name of input file
        FILE=${filename}
        echo ${FILE}
        # # Path where the file is located
        # INPUT_FILE=${PATH_FILE}${FILE}
        # # Each line in the file is an argument passed to another file.
        # # Here the lines are extracted and stored as argument/variable.
        # mkdir_path=$(sed '1q;d' ${INPUT_FILE})
        # arg_mutect2=$(sed '2q;d' ${INPUT_FILE})
        # file_output=$(sed '3q;d' ${INPUT_FILE})
        # # The file from which a function is called.
        # source /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/make_vcf_auto.sh
        # # The function called from the file.
        # # The three arguments are given.
        # # The argument in quotes "" consists of multiple words (split by space), 
        # # so it must be passed in quotes. Otherwise, only the first 'word' 
        # # will be grabbed before the first space character.
        # mutect2_vcf ${mkdir_path} "${arg_mutect2}" ${file_output}
        echo "EIND"
    done    
done

