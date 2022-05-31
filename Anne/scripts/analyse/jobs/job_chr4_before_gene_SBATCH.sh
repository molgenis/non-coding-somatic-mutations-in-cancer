#!/usr/bin/bash
#SBATCH --job-name='job_chr4_before_gene'
#SBATCH --output='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/analyse/jobs/job_chr4_before_gene.out'
#SBATCH --error='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/analyse/jobs/job_chr4_before_gene.err'
#SBATCH --time='100:59:59'
#SBATCH --cpus-per-task='8'
#SBATCH --mem='64'GB
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp='10gb'
ml Miniconda3
source activate stage2
python3 '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/analyse/before_gene_multi.py' 'chr4'
