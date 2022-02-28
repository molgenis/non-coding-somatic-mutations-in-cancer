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

# Path to directory with files
# PATH_DIR='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/merge_H/'
PATH_DIR='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/merge_B/'
# Path to script
PATH_SCRIPT='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/new_plan/'

FILE_LIST=${PATH_DIR}file_list.txt
rm -f ${FILE_LIST}

ml BCFtools/1.11-GCCcore-7.3.0

# Loop over files and get the file who ends with .vcf.gz
for filename in ${PATH_DIR}*.vcf.gz; do
    echo '-------------------------------------'
    echo "$filename"
    # Get basename
    BASENAME=$( echo "$(basename -- $filename)")
    # Get code (before .)
    IFS='.' read -r codeName string <<< "$BASENAME"
    # Rename column names
    python3 ${PATH_SCRIPT}change_colname_newdata.py ${filename} ${PATH_DIR}
    # Check if nohead.vcf exists (made by python file)
    if [ -f "${PATH_DIR}nohead.vcf" ]; then
        # Get header
        bcftools view --header-only $filename | sed 's/##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/' | sed 's/##contig=<ID=MT/##contig=<ID=chrM/' > ${PATH_DIR}header.vcf
        # Remove last sentence
        sed -i '$ d' ${PATH_DIR}header.vcf
        #no headers
        #bcftools view --header-only $filename | sed 's/##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/' | sed 's/##contig=<ID=MT/##contig=<ID=chrM/' > ${PATH_DIR}nohead.vcf 
        #header en no header merge
        cat ${PATH_DIR}header.vcf ${PATH_DIR}nohead.vcf >> ${PATH_DIR}merge_${codeName}.vcf
        #https://github.com/samtools/bcftools/issues/668
        bcftools view ${PATH_DIR}merge_${codeName}.vcf -Oz -o ${PATH_DIR}merge_${codeName}.vcf.gz
        bcftools index ${PATH_DIR}merge_${codeName}.vcf.gz
        # Make file list for bcftools merge
        echo "${PATH_DIR}merge_${codeName}.vcf.gz" >> "${FILE_LIST}"
        # Remove files
        rm ${PATH_DIR}header.vcf
        rm ${PATH_DIR}nohead.vcf
        rm ${PATH_DIR}merge_${codeName}.vcf
    fi
       
done

ml BCFtools/1.11-GCCcore-7.3.0
#tr '\n' ' ' < ${FILE_LIST}
#MERGE_COM=$(sed '1q;d' ${FILE_LIST})
#echo ${MERGE_COM}
# Merge vcf files
bcftools merge  --file-list ${FILE_LIST}  -O z -o ${PATH_DIR}mergAll_merged.vcf.gz
bcftools index ${PATH_DIR}mergAll_merged.vcf.gz
# Remove files
rm ${PATH_DIR}merge_*.vcf.gz
rm ${PATH_DIR}merge_*.vcf.gz.csi