#!/usr/bin/bash

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

source /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/read_yaml.sh
YAML_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/config.yaml

SCRIPT_PATH=$(yaml $YAML_PATH "['database_scripts']")
PATH_DATA=$(yaml $YAML_PATH "['database_path']")
# DB_NAME='00Database_internship_UPDATE2.0.db'
DATABASE_SNP=$(yaml $YAML_PATH "['database']") #${DB_NAME}.db
GENE_FILE=$(yaml $YAML_PATH "['checkGene']")

#A, D, H, W, M, U, E, G, K, N, O, T, R, S, C, B, P, L
# GEDAAN: 'O' 'A' 'D' 'H' 'W' 'U' 'G' 'K' 'N' 'T' 'R' 'B' 'C'
LETTERS1=('LIRI-JP' 'LMS-FR' 'LUAD-US')
# LETTERS1=('A' 'D' 'H' 'W' 'U' 'G' 'K' 'N' 'O' 'T' 'R')
# LETTERS2=('B' 'C')
# LETTERS3=('L')
# LETTERS4=('P')
# LETTERS5=('S')
# LETTERS6=('E')
# LETTERS7=('M')

python3 ${SCRIPT_PATH}create_db.py ${DATABASE_SNP}

for LETTER in "${LETTERS1[@]}"
do
    for filename in ${PATH_DATA}${LETTER}*_db_NEW.tsv.gz
    do
        # DB_NAME="${filename}
        # DATABASE_SNP=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/"$(b=${filename##*/}; echo ${b%%.*})".db
        # python3 ${SCRIPT_PATH}create_db.py ${DATABASE_SNP}
        # Name of input file
        echo "$(b=${filename##*/}; echo ${b%%.*})"
        echo '------------ fill DB'
        python3 ${SCRIPT_PATH}fill_db.py ${DATABASE_SNP} ${filename}
        # echo '------------ Check genes'
        # python3 ${SCRIPT_PATH}check_gene.py ${DATABASE_GENE} ${GENE_FILE}
        echo "EIND" "$(b=${filename##*/}; echo ${b%%.*})"
    done
done    

echo 'END END END'


# python3 ${PATH_DB}check_gene_filedb.py