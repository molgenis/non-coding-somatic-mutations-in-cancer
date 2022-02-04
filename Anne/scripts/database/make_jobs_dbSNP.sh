#!/usr/bin/bash

#SBATCH --job-name=db_test
#SBATCH --output=db_test.out
#SBATCH --error=db_test.err
#SBATCH --time=159:59:59
#SBATCH --cpus-per-task=20
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

# Make array of chromosomes
#seq FIRST STEP LAST
chrom_num=($(seq 1 1 22))
chrom_num+=("X" "Y")

# The name of the database
DB_NAME='Database_internship_gene_long_NEW2.0'
# Path where the database is/should be stored
DB_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/${DB_NAME}.db
# Path where the files should be stored
DB_FILES=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/files_to_annotate/
# Path to the scripts
SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/all_git/non-coding-somatic-mutations-in-cancer/Anne/scripts/database/
CPUS=20


for i in "${chrom_num[@]}"
do
    # Path to the scripts
    JOB_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/all_git/non-coding-somatic-mutations-in-cancer/Anne/scripts/database/jobs/
    # File
    DB_SNP_FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/files_to_annotate/annotate/chr${i}_ann.vcf
    FILE=${JOB_PATH}chr${i}_create_vcf_and_annotate_dbSNP.sh
    echo $i
    echo "#!/usr/bin/bash" > ${FILE}
    echo "" >> ${FILE}
    echo "#SBATCH --job-name=file${i}" >> ${FILE}
    echo "#SBATCH --output=${JOB_PATH}file${i}.out" >> ${FILE}
    echo "#SBATCH --error=${JOB_PATH}file${i}.err" >> ${FILE}
    echo "#SBATCH --time=159:59:59" >> ${FILE}
    echo "#SBATCH --cpus-per-task=${CPUS}" >> ${FILE}
    echo "#SBATCH --mem=96gb" >> ${FILE}
    echo "#SBATCH --nodes=1" >> ${FILE}
    echo "#SBATCH --open-mode=append" >> ${FILE}
    echo "#SBATCH --export=NONE" >> ${FILE}
    echo "#SBATCH --get-user-env=L" >> ${FILE}
    echo "" >> ${FILE}
    echo "# Load ml" >> ${FILE}
    echo "ml BCFtools/1.11-GCCcore-7.3.0" >> ${FILE}
    echo "ml Anaconda3/5.3.0" >> ${FILE}
    echo "source activate stage" >> ${FILE}
    echo "" >> ${FILE}
    echo "python3 ${SCRIPT_PATH}create_vcf_file.py ${DB_PATH} ${DB_FILES} ${i}" >> ${FILE}
    echo "bgzip ${DB_FILES}chr${i}_db.vcf #bcftools view file.vcf -Oz -o file.vcf.gz" >> ${FILE}
    echo "tabix ${DB_FILES}chr${i}_db.vcf.gz #bcftools index file.vcf.gz" >> ${FILE}
    echo "bcftools annotate -c ID -a /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/dbSNP/per_chr/chr${i}_merge_All_20180423.vcf.gz   -o ${DB_FILES}annotate/chr${i}_ann.vcf  ${DB_FILES}chr${i}_db.vcf.gz" >> ${FILE}
    if  [[ ${i} == 1 ]]
    then
        echo "------${i}"
        echo "python3 ${SCRIPT_PATH}filter_file_dbSNP.py ${DB_PATH} ${DB_SNP_FILE} ALTER" >> ${FILE}
    else
        echo "python3 ${SCRIPT_PATH}filter_file_dbSNP.py ${DB_PATH} ${DB_SNP_FILE} EMPTY" >> ${FILE}
    fi    
    
    echo "" >> ${FILE}
    echo "RUN ${FILE}"
    source ${FILE}
done

