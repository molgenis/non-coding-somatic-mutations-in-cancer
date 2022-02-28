#!/usr/bin/bash

#SBATCH --job-name=cha_col
#SBATCH --output=cha_col.out
#SBATCH --error=cha_col.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Anaconda3/5.3.0
source activate stage

# PATH_DIR='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/merge_B/'
PATH_DIR='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/merge_B/'
PATH_NEW='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/merge_B/'
PATH_SCRIPT='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/new_plan/'
rm -f ${PATH_SCRIPT}file.list

for filename in ${PATH_DIR}*.vcf.gz; do
    echo "$filename"
    # Get basename
    BASENAME=$( echo "$(basename -- $filename)")
    # Get code (before .)
    IFS='.' read -r codeName string <<< "$BASENAME"
    # Get header: https://www.biostars.org/p/352580/
    cat $filename | grep '^##' | sed 's/=/\,/g' | sed 's/#//g' > ${PATH_NEW}header.txt
    python3 ${PATH_SCRIPT}change_colname_newdata.py ${filename} ${PATH_NEW}
    # header en no header plakken
    cat ${PATH_NEW}header.txt ${PATH_NEW}nohead.vcf >> ${PATH_NEW}merge.vcf
    gzip ${PATH_NEW}merge_${codeName}.vcf
    echo "${PATH_NEW}merge_${codeName}.vcf.gz" >> ${PATH_SCRIPT}file.list

    rm ${PATH_NEW}header.txt
    rm ${PATH_NEW}nohead.vcf
    rm ${PATH_NEW}merge_${codeName}.vcf
done

ml BCFtools/1.11-GCCcore-7.3.0
bcftools merge  --file-list file.list  -O z -o ${PATH_NEW}All_merged.vcf.gz
# rm ${PATH_NEW}merge_*.vcf.gz