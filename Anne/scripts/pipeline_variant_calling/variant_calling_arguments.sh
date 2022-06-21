#!/usr/bin/bash

# Gets the arguments for running Mutect from files and forwards them to 
# another bash file that Mutect2 eventually runs

#SBATCH --job-name=calling
#SBATCH --output=calling.out
#SBATCH --error=calling.err
#SBATCH --time=149:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

echo "job vcf ${METHOD}"

# The path where the file is located
PATH_FILE=${GENERAL_PATH}${SAMPLE}/${CHROM}/${TYPE_ALN}/
# Runs over all files matching the path and ending in .txt
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
    source ${SCRIPT_PATH}variant_calling.sh
    # The function called from the file.
    # The three arguments are given.
    # The argument in quotes "" consists of multiple words (split by space), 
    # so it must be passed in quotes. Otherwise, only the first 'word' 
    # will be grabbed before the first space character.
    mutect2_vcf ${mkdir_path} "${arg_mutect2}" ${file_output}
    echo "EIND job vcf ${METHOD} - ${filename}"
done
echo "EIND job vcf ${METHOD} - ${SAMPLE}"    

