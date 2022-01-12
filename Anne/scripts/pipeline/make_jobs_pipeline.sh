#!/usr/bin/bash

# The samples
SAMPLES=(S1 S2 S3 S4 S5 S6)
# The tissues of the samples
TISSUE=( '4094  4099  4104' '4109  4113  4114  4118  4119' '4123  4124  4129' '4128  4133  4134' '4139  4138  5041' '5043  5042  5044')

# The choosen align method
METHOD_BIG=bwa_mem #bowtie,   bwa_aln, bwa_mem
CHROM_CHOOSEN=chr21

# Passes different things to the variable VARIABLE and METH_FILE_BIG, 
# depending on the variable METHOD_BIG
# VARIABLE: Scripts names needed
if [ "${METHOD_BIG}" == "bwa_aln" ]; then
    # VARIABLE=( job_align_aln job_vcf_aln )
    METH_FILE_BIG=aln #bowtie2, aln, mem
elif [ "${METHOD_BIG}" == "bwa_mem" ]; then
    # VARIABLE=( job_align_mem job_vcf_mem ) 
    METH_FILE_BIG=mem #bowtie2, aln, mem   
elif [ "${METHOD_BIG}" == "bowtie" ]; then
    # VARIABLE=( job_align_bowtie job_vcf_bowtie )
    METH_FILE_BIG=bowtie2 #bowtie2, aln, mem  
else
    echo "ERROR"
fi

# Make new directory
mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/pipeline/jobs/

# Loop over samples
for i in "${!SAMPLES[@]}"
do  
    # Path to the jobs
    JOB_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/pipeline/jobs/
    # Filename of the new jobs
    FILE_NAME_JOB=${JOB_PATH}"${SAMPLES[i]}"_${METHOD_BIG}
    JOB_FILE=${FILE_NAME_JOB}.sh
    echo "${TISSUE[i]}"
    # Make array of variable TISSUE
    IFS=', ' read -r -a array <<< "${TISSUE[i]}"
    # Everything that comes after this will be in the new job
    echo "#!/usr/bin/bash" > "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo "#SBATCH --job-name=${FILE_NAME_JOB}" >> "${JOB_FILE}"
    echo "#SBATCH --output=${FILE_NAME_JOB}.out" >> "${JOB_FILE}"
    echo "#SBATCH --error=${FILE_NAME_JOB}.err" >> "${JOB_FILE}"
    echo "#SBATCH --time=89:59:59" >> "${JOB_FILE}"
    echo "#SBATCH --cpus-per-task=8" >> "${JOB_FILE}"
    echo "#SBATCH --mem=96gb" >> "${JOB_FILE}"
    echo "#SBATCH --nodes=1" >> "${JOB_FILE}"
    echo "#SBATCH --open-mode=append" >> "${JOB_FILE}"
    echo "#SBATCH --export=NONE" >> "${JOB_FILE}"
    echo "#SBATCH --get-user-env=L" >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo "# Load packages" >> "${JOB_FILE}"
    echo "ml BCFtools/1.11-GCCcore-7.3.0" >> "${JOB_FILE}"
    echo "ml Bowtie2/2.3.4.2-foss-2018b" >> "${JOB_FILE}"
    echo "ml BWA/0.7.17-GCCcore-7.3.0" >> "${JOB_FILE}"
    echo "ml cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare" >> "${JOB_FILE}"
    echo "ml FastQC/0.11.8-Java-11-LTS" >> "${JOB_FILE}"
    echo "ml GATK/4.1.4.1-Java-8-LTS" >> "${JOB_FILE}"
    echo "ml Java/11-LTS" >> "${JOB_FILE}"
    echo "#ml libjpeg-turbo/2.0.2-GCCcore-7.3.0 #?" >> "${JOB_FILE}"
    echo "ml picard/2.20.5-Java-11-LTS" >> "${JOB_FILE}"
    #
    echo "" >> "${JOB_FILE}"
    echo "#Chrom" >> "${JOB_FILE}"
    echo "CHROM=${CHROM_CHOOSEN}" >> "${JOB_FILE}"
    echo "# The path to the file of the reference genome that will be used." >> "${JOB_FILE}"
    echo 'if [ "${CHROM}" != "chrall" ]; then' >> "${JOB_FILE}"
    echo "  PATH_GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/" >> "${JOB_FILE}"
    echo "else" >> "${JOB_FILE}"
    echo "  PATH_GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/ #chrall" >> "${JOB_FILE}"
    echo "fi" >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo 'GENOOM=${PATH_GENOOM}${CHROM}.fa' >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo "ml SAMtools/1.9-foss-2018b" >> "${JOB_FILE}"
    echo "# When file doesn't exist." >> "${JOB_FILE}"
    echo 'if [ ! -f ${PATH_GENOOM}${CHROM}.dict ]; then' >> "${JOB_FILE}"
    echo "  echo 'FILE does not exist.'" >> "${JOB_FILE}"
    echo "  # dict file" >> "${JOB_FILE}"
    echo '  gatk CreateSequenceDictionary -R ${GENOOM}' >> "${JOB_FILE}"
    echo "  # index file" >> "${JOB_FILE}"
    echo '  samtools faidx ${GENOOM}' >> "${JOB_FILE}"
    echo "fi" >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    # INSTALL PARAMS
    echo "# Array of tissue numbers" >> "${JOB_FILE}"
    echo "TISSUE_ARR=("${array[@]}")" >> "${JOB_FILE}"
    echo "# Sample numbers" >> "${JOB_FILE}"
    echo "SAMPLE="${SAMPLES[i]}"" >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo "METHOD=${METHOD_BIG} #bowtie,   bwa_aln, bwa_mem" >> "${JOB_FILE}"
    echo "METH_FILE=${METH_FILE_BIG} #bowtie2, aln, mem" >> "${JOB_FILE}"
    echo 'TYPE_ALN=mutect_${METHOD}' >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo "GENERAL_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/" >> "${JOB_FILE}"
    echo "SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/pipeline/" >> "${JOB_FILE}"
    echo "# Path to the file of the panel of normals (PoN)" >> "${JOB_FILE}"
    echo "PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/merge_somatic-b37_Mutect2-WGS-panel-b37.vcf" >> "${JOB_FILE}"
    echo "# Path to the file of the germline resource (GR)" >> "${JOB_FILE}"
    echo "GR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GermlineResource/merge_af-only-gnomad.raw.sites.vcf" >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo "# for python file parameters: automatic_script_ob.py" >> "${JOB_FILE}"
    echo "# Number of tumors you want to combine while running Mutect2." >> "${JOB_FILE}"
    echo "NUMBER_OF_TUMORS_py=1" >> "${JOB_FILE}"
    echo "# Number of hc you want to combine while running Mutect2" >> "${JOB_FILE}"
    echo "NUMBER_OF_HC_py=1" >> "${JOB_FILE}"
    echo "# What type of tumor you want to have ('both', 'tFL' or 'FL')" >> "${JOB_FILE}"
    echo "TYPE_SAMPLE_py='both'" >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    # PIPELINE
    echo "# load conda and activate to execute python script" >> "${JOB_FILE}"
    echo "ml Anaconda3/5.3.0" >> "${JOB_FILE}"
    echo "source activate stage" >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo "# RUN: file_prep.sh" >> "${JOB_FILE}"
    echo '#source ${SCRIPT_PATH}file_prep.sh' >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo 'echo ${METHOD}' >> "${JOB_FILE}"
    echo "# RUN: align.sh" >> "${JOB_FILE}"
    echo "echo 'align.sh'" >> "${JOB_FILE}"
    echo 'source ${SCRIPT_PATH}align.sh' >> "${JOB_FILE}"
    echo 'mkdir -p ${GENERAL_PATH}${SAMPLE}/${CHROM}/mutect_${METHOD}/' >> "${JOB_FILE}"
    echo "# RUN: make_file_vc.py" >> "${JOB_FILE}"
    echo "echo 'make_file_vc.py'" >> "${JOB_FILE}"
    echo 'python3 ${SCRIPT_PATH}make_file_vc.py ${GENERAL_PATH} ${NUMBER_OF_TUMORS_py} ${NUMBER_OF_HC_py} ${TYPE_SAMPLE_py} ${METHOD} ${METH_FILE} ${CHROM}' >> "${JOB_FILE}"
    echo "# RUN: change_sample_name.sh" >> "${JOB_FILE}"
    echo "echo 'change_sample_name'" >> "${JOB_FILE}"
    echo 'source ${SCRIPT_PATH}change_sample_name.sh' >> "${JOB_FILE}"
    echo "# RUN: variant_calling_arguments.sh" >> "${JOB_FILE}"
    echo "# mutect2" >> "${JOB_FILE}"
    echo "echo 'variant_calling_arguments.sh'" >> "${JOB_FILE}"
    echo 'source ${SCRIPT_PATH}variant_calling_arguments.sh' >> "${JOB_FILE}"
    echo "# RUN: compare_vcf.sh" >> "${JOB_FILE}"
    echo "echo 'compare_vcf.sh'" >> "${JOB_FILE}"
    echo 'source ${SCRIPT_PATH}compare_vcf.sh' >> "${JOB_FILE}"
    echo "# RUN: merge_vcf.sh" >> "${JOB_FILE}"
    echo "# Merge files" >> "${JOB_FILE}"
    echo "COMP_TYPE=mutect2 # manual, mutect2" >> "${JOB_FILE}"
    echo "echo 'merge_vcf1'" >> "${JOB_FILE}"
    echo 'source ${SCRIPT_PATH}merge_vcf.sh' >> "${JOB_FILE}"
    echo "COMP_TYPE=manual # manual, mutect2" >> "${JOB_FILE}"
    echo "echo 'vmerge_vcf2'" >> "${JOB_FILE}"
    echo 'source ${SCRIPT_PATH}merge_vcf.sh' >> "${JOB_FILE}"
    echo "# RUN: dbSNP_annotate.sh" >> "${JOB_FILE}"
    echo "# filteren van dbSNP (annoteren)" >> "${JOB_FILE}"
    echo "echo 'dbSNP_annotate'" >> "${JOB_FILE}"
    echo 'source ${SCRIPT_PATH}dbSNP_annotate.sh' >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo "" >> "${JOB_FILE}"
    echo "echo 'THE END END END END'" >> "${JOB_FILE}"
done