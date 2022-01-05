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

#seq FIRST STEP LAST
chrom_num=($(seq 1 1 22))
chrom_num+=("X" "Y")

DB_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene.db
DB_FILES=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/files_to_annotate/
SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/all_git/non-coding-somatic-mutations-in-cancer/Anne/scripts/database/

for i in "${chrom_num[@]}"
do
    JOB_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/all_git/non-coding-somatic-mutations-in-cancer/Anne/scripts/database/jobs/
    FILE=${JOB_PATH}file${i}.sh
    echo $i
    echo "#!/usr/bin/bash" > ${FILE}
    echo "" >> ${FILE}
    echo "#SBATCH --job-name=file${i}" >> ${FILE}
    echo "#SBATCH --output=${JOB_PATH}file${i}.out" >> ${FILE}
    echo "#SBATCH --error=${JOB_PATH}file${i}.err" >> ${FILE}
    echo "#SBATCH --time=159:59:59" >> ${FILE}
    echo "#SBATCH --mem=96gb" >> ${FILE}
    echo "#SBATCH --nodes=1" >> ${FILE}
    echo "#SBATCH --open-mode=append" >> ${FILE}
    echo "#SBATCH --export=NONE" >> ${FILE}
    echo "#SBATCH --get-user-env=L" >> ${FILE}
    echo "" >> ${FILE}
    echo "# Load ml" >> ${FILE}
    echo "ml Anaconda3/5.3.0" >> ${FILE}
    echo "source activate stage" >> ${FILE}
    echo "" >> ${FILE}
    echo "python3 ${SCRIPT_PATH}create_vcf_file.py ${DB_PATH} ${DB_FILES} ${i}" >> ${FILE}
    echo "bgzip ${DB_FILES}chr${i}_db.vcf #bcftools view file.vcf -Oz -o file.vcf.gz" >> ${FILE}
    echo "tabix ${DB_FILES}chr${i}_db.vcf.gz #bcftools index file.vcf.gz" >> ${FILE}
    echo "" >> ${FILE}
    echo "" >> ${FILE}
done