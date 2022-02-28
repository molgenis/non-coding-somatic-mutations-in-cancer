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

# PATH_DIR='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/merge_H/'
PATH_DIR='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/merge_B/'
PATH_NEW='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/merge_B/'
PATH_SCRIPT='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/new_plan/'

FILE_LIST=${PATH_DIR}file_list.txt
rm -f ${FILE_LIST}

#echo "" > "${FILE_LIST}"

for filename in ${PATH_DIR}0*.vcf.gz; do
    echo '-------------------------------------'
    echo "$filename"
    # Get basename
    BASENAME=$( echo "$(basename -- $filename)")
    # Get code (before .)
    IFS='.' read -r codeName string <<< "$BASENAME"
    python3 ${PATH_SCRIPT}change_colname_newdata.py ${filename} ${PATH_DIR}
    
    if [ -f "${PATH_DIR}nohead.vcf" ]; then
        # Get header: https://www.biostars.org/p/352580/
        cat $filename | grep '^##' | sed 's/=/\,/g' | sed 's/#//g' > ${PATH_DIR}header.txt
    
        # header en no header plakken
        cat ${PATH_DIR}header.txt ${PATH_DIR}nohead.vcf >> ${PATH_DIR}merge_${codeName}.vcf
        gzip ${PATH_DIR}merge_${codeName}.vcf
        echo "${PATH_DIR}merge_${codeName}.vcf.gz" >> "${FILE_LIST}"
        
        rm ${PATH_DIR}header.txt
        rm ${PATH_DIR}nohead.vcf
        #rm ${PATH_DIR}merge_${codeName}.vcf
    fi
       
done

ml BCFtools/1.11-GCCcore-7.3.0
tr '\n' ' ' < ${FILE_LIST}
MERGE_COM=$(sed '1q;d' ${FILE_LIST})
echo ${MERGE_COM}
#bcftools merge  --file-list ${FILE_LIST}  -O z -o ${PATH_DIR}All_merged.vcf.gz
bcftools merge --merge all ${MERGE_COM} -O z -o ${PATH_DIR}mergAll_merged.vcf.gz --force-samples
# rm ${PATH_DIR}merge_*.vcf.gz