#!/usr/bin/bash

#SBATCH --job-name=gene_db
#SBATCH --output=gene_db.out
#SBATCH --error=gene_db.err
#SBATCH --time=89:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage

PATH_DB=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/all_git/non-coding-somatic-mutations-in-cancer/Anne/scripts/database/

python3 ${PATH_DB}check_gene_filedb.py