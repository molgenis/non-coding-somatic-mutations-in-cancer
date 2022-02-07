#!/usr/bin/bash

#SBATCH --job-name=db_test
#SBATCH --output=db_test.out
#SBATCH --error=db_test.err
#SBATCH --time=99:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage

SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/database2/
PATH_DATA=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/data_db/ #todo /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_copy/
DB_NAME='Database_internship_gene_long_NEW_NEW_UPDATE'
DATABASE_GENE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/${DB_NAME}.db
GENE_FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/snp132_ucsc_hg19_checkGene.bed

#A, D, H, W, M, U, E, G, K, N, O, T, R, S, C, B, P, L
LETTERS1=('ALL')
# LETTERS1=('A' 'D' 'H' 'W' 'U' 'G' 'K' 'N' 'O' 'T' 'R')
# LETTERS2=('B' 'C')
# LETTERS3=('L')
# LETTERS4=('P')
# LETTERS5=('S')
# LETTERS6=('E')
# LETTERS7=('M')

python3 ${SCRIPT_PATH}create_db.py ${DATABASE_GENE}

for LETTER in "${LETTERS1[@]}"
do
    for filename in ${PATH_DATA}${LETTER}*dbNEW.tsv
    do
        # Name of input file
        echo ${filename}
        echo '------------ fill DB'
        python3 ${SCRIPT_PATH}fill_db.py ${DATABASE_GENE} ${filename}
        # echo '------------ Check genes'
        # python3 ${SCRIPT_PATH}check_gene.py ${DATABASE_GENE} ${GENE_FILE}
        echo "EIND" ${filename}
    done
done    

echo 'END END END'


# python3 ${PATH_DB}check_gene_filedb.py