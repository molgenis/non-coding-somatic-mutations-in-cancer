#!/usr/bin/bash

#SBATCH --job-name=db_test
#SBATCH --output=db_test.out
#SBATCH --error=db_test.err
#SBATCH --time=159:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

# Load modules
ml Anaconda3/5.3.0
source activate stage

# Path to the scripts
SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/database/
# Path to the data to be read
PATH_DATA=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/data_db/
# The name of the database
DB_NAME='Database_internship_gene_long_NEW2.0'
# Path where the database is/should be stored
DATABASE_GENE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/${DB_NAME}.db
# File with the genes (chromosome - start - end - etc)
GENE_FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/snp132_ucsc_hg19_checkGene.bed

# Array with letters
LETTERS1=('A' 'D' 'H' 'W' 'U' 'G' 'K' 'N' 'O' 'T' 'R')
LETTERS2=('B' 'C' 'L' 'P' 'S' 'E')
LETTERS3=('M')

# Run create_db.py
# python3 ${SCRIPT_PATH}create_db.py ${DATABASE_GENE}

for LETTER in "${LETTERS1[@]}"
do
    # Loop over all files that ends with db.tsv
    for filename in ${PATH_DATA}${LETTER}*db.tsv
    do
        # Name of input file
        echo ${filename}
        echo '------------ fill DB'
        # Run fill_db.py
        python3 ${SCRIPT_PATH}fill_db.py ${DATABASE_GENE} ${filename}
        # echo '------------ Check genes'
        # Run check_gene.py
        # python3 ${SCRIPT_PATH}check_gene.py ${DATABASE_GENE} ${GENE_FILE}
        echo "EIND" ${filename}
    done
done    

echo 'END END END'


# python3 ${PATH_DB}check_gene_filedb.py