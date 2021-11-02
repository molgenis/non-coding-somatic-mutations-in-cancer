#!/usr/bin/bash
# Name of input file
FILE=bowtie_S1_numT_1_numHC_1_SS6004094_SS6004099.txt
# Path where the file is located
INPUT_FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/S1/${FILE}
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
