#!/usr/bin/bash

# Walks over the projects and calls a python script that populates the database.

#SBATCH --job-name=db_test
#SBATCH --output=db_test.out
#SBATCH --error=db_test.err
#SBATCH --time=120:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage

#Path to scripts
SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/database/
#Path to database forder
PATH_DATA=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_data_db/
#Path to database
DATABASE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/db_laatste_copy.db

#A, D, H, W, M, U, E, G, K, N, O, T, R, S, C, B, P, L
# GEDAAN: 'O' 'A' 'D' 'H' 'W' 'U' 'G' 'K' 'N' 'T' 'R' 'B' 'C'
LETTERS=('A' 'D' 'H' 'W' 'U' 'G' 'K' 'N' 'O' 'T' 'R' 'B' 'C' 'L' 'P' 'S' 'E' 'M')

python3 ${SCRIPT_PATH}create_db.py ${DATABASE}

for LETTER in "${LETTERS[@]}"
do
    for filename in ${PATH_DATA}${LETTER}*_db_NEW.tsv.gz
    do
        # RUN: create_db.py
        python3 ${SCRIPT_PATH}create_db.py ${DATABASE}
        # Name of input file
        echo "$(b=${filename##*/}; echo ${b%%.*})"
        echo '------------ fill DB'
        # RUN: fill_db.py
        python3 ${SCRIPT_PATH}fill_db.py ${DATABASE} ${filename}
        echo "EIND" "$(b=${filename##*/}; echo ${b%%.*})"
    done
done    

echo 'END END END'