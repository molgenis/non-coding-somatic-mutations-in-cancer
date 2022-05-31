TYPE_DATA=('before_gene' 'DNase' 'TFBS' 'UCNE' \
    'genotype' 'per_snp' 'region')

for type_data in "${TYPE_DATA[*]}"
do
    PYTHON_SCRIPT='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/analyse/'${type_data}'_multi.py'
    #PROCESS_LANE_SCRIPT_LOC='/groups/umcg-weersma/tmp04/projects/lpmc_v2/ongoing/seurat_preprocess_samples/scripts/lpmcv2_lane_to_seurat.R'
    JOB_DIR='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/analyse/jobs'
    #JOB_DIR='/groups/umcg-weersma/tmp04/projects/lpmc_v2/ongoing/seurat_preprocess_samples/jobs/'
    # these are the lanes to run through cellranger
    CHROMS=('chr1' 'chr2' 'chr3' 'chr4' \
    'chr5' 'chr6' 'chr7' 'chr8' \
    'chr9' 'chr10' 'chr11' 'chr12' \
    'chr13' 'chr14' 'chr15' 'chr16' \
    'chr17' 'chr18' 'chr19' 'chr20' \
    'chr21' 'chr22' 'chrX' 'chrY' )
    CORES='8'
    MEMORY_GB='64'
    TMP_SIZE='10gb'
    RUNTIME='100:59:59'


    # check each run
    for chr in "${CHROMS[*]}"
    do
        JOB_NAME='read_'${chr}'_'${type_data}
        JOB_LOC=${JOB_DIR}'/'${JOB_NAME}'_SBATCH.sh'
        JOB_OUT=${JOB_DIR}'/'${JOB_NAME}'.out'
        JOB_ERR=${JOB_DIR}'/'${JOB_NAME}'.err'
            # echo the header
            echo '#!/bin/bash
        #SBATCH --job-name='${JOB_NAME}'
        #SBATCH --output='${JOB_OUT}'
        #SBATCH --error='${JOB_ERR}'
        #SBATCH --time='${RUNTIME}'
        #SBATCH --cpus-per-task='${CORES}'
        #SBATCH --mem='${MEMORY_GB}'GB
        #SBATCH --nodes=1
        #SBATCH --export=NONE
        #SBATCH --get-user-env=L
        #SBATCH --tmp='${TMP_SIZE}'
        ml Miniconda3
        source activate stage2
        python3 '${PYTHON_SCRIPT}' '${chr}' ' >> ${JOB_LOC}
    done
done