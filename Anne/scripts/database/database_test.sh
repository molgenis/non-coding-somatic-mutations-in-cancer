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

ml Anaconda3/5.3.0
source activate stage

PATH_DB=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/all_git/non-coding-somatic-mutations-in-cancer/Anne/scripts/database/
PATH_DATA=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/data_db/


for filename in ${PATH_DATA}B*.tsv
do
    # Name of input file
    echo ${filename}
    python3 ${PATH_DB}database_test.py ${filename} 'db_test'
    echo "EIND" ${filename}
done    

echo 'END END END'


# python3 ${PATH_DB}check_gene_filedb.py