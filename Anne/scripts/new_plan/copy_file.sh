#!/usr/bin/bash

#SBATCH --job-name=AF
#SBATCH --output=AF.out
#SBATCH --error=AF.err
#SBATCH --time=80:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=116gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage

#PATH_OLD=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/B/
#PATH_NEW=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/All_B/

PATH_OLD=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/H/
PATH_NEW=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/All_H/


for FOL in ${PATH_OLD}*/ ; do
    echo "$FOL"
    echo 'HOI'
    for filename in ${FOL}*.vcf.gz*; do
        echo "$filename"
        BASENAME=$( echo "$(basename -- $filename)")
        echo ${BASENAME}
        cp ${filename} ${PATH_NEW}${BASENAME}
    done
done