TYPE_DATA=(before_gene DNase TFBS UCNE genotype per_snp region)

for type_data in "${!TYPE_DATA[@]}"
do
    echo $type_data
    PYTHON_SCRIPT='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/analyse_layers/'${TYPE_DATA[type_data]}'_multi.py'
    #PROCESS_LANE_SCRIPT_LOC='/groups/umcg-weersma/tmp04/projects/lpmc_v2/ongoing/seurat_preprocess_samples/scripts/lpmcv2_lane_to_seurat.R'
    JOB_DIR='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/analyse_layers/jobs'
    #JOB_DIR='/groups/umcg-weersma/tmp04/projects/lpmc_v2/ongoing/seurat_preprocess_samples/jobs/'
    ANALYSE='/groups/umcg-wijmenga/tmp04/projects/lude_vici_2021/rawdata/cancer_data/analyse_layers/'
    # these are the lanes to run through cellranger
    CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
    CORES='8'
    MEMORY_GB='50'
    TMP_SIZE='10gb'
    RUNTIME='160:59:59'


    # check each run
    for chr in "${!CHROMS[@]}"
    do
        echo $chr
        JOB_NAME='job_'${CHROMS[chr]}'_'${TYPE_DATA[type_data]}
        JOB_LOC=${JOB_DIR}'/'${JOB_NAME}'_SBATCH.sh'
        JOB_OUT=${JOB_DIR}'/'${JOB_NAME}'.out'
        JOB_ERR=${JOB_DIR}'/'${JOB_NAME}'.err'
        # echo the header
        echo "#!/usr/bin/bash" > "${JOB_LOC}"
        echo "#SBATCH --job-name='${JOB_NAME}'" >> "${JOB_LOC}"
        echo "#SBATCH --output='${JOB_OUT}'" >> "${JOB_LOC}"
        echo "#SBATCH --error='${JOB_ERR}'" >> "${JOB_LOC}"
        echo "#SBATCH --time='${RUNTIME}'" >> "${JOB_LOC}"
        echo "#SBATCH --cpus-per-task='${CORES}'" >> "${JOB_LOC}"
        echo "#SBATCH --mem='${MEMORY_GB}'GB" >> "${JOB_LOC}"
        echo "#SBATCH --nodes=1" >> "${JOB_LOC}"
        echo "#SBATCH --export=NONE" >> "${JOB_LOC}"
        echo "#SBATCH --get-user-env=L" >> "${JOB_LOC}"
        echo "#SBATCH --tmp='${TMP_SIZE}'" >> "${JOB_LOC}"
        # echo "cd ${TMPDIR}" >> "${JOB_LOC}"
        # echo "NEW_PATH='${CHROMS[chr]}'_'${TYPE_DATA[type_data]}'" >> "${JOB_LOC}"
        # echo "mkdir -p ${NEW_PATH}" >> "${JOB_LOC}"
        echo "ml Miniconda3" >> "${JOB_LOC}"
        echo "source activate stage2" >> "${JOB_LOC}"        
        echo "python3 '${PYTHON_SCRIPT}' '${CHROMS[chr]}'" >> "${JOB_LOC}" #${TMPDIR}'${NEW_PATH}'
        # echo "cp -r '${NEW_PATH}'*.tsv" >> "${JOB_LOC}"
        # cp -r ${CHROMS[chr]}_${TYPE_DATA[type_data]} 
    done
done