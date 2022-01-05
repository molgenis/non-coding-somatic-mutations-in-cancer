#!/usr/bin/bash

SAMPLES=(S1 S2 S3 S4 S5 S6)
TISSUE=( '4094  4099  4104' '4109  4113  4114  4118  4119' '4123  4124  4129' '4128  4133  4134' '4139  4138  5041' '5043  5042  5044')

METHOD_BIG=bwa_aln #bowtie,   bwa_aln, bwa_mem
CHROM_CHOOSEN=chr21

if [ "${METHOD}" == "bwa_aln" ]; then
    VARIABLE=( job_align_aln job_vcf_aln )
    METH_FILE_BIG=aln #bowtie2, aln, mem
elif [ "${METHOD}" == "bwa_mem" ]; then
    VARIABLE=( job_align_mem job_vcf_mem ) 
    METH_FILE_BIG=mem #bowtie2, aln, mem   
elif [ "${METHOD}" == "bowtie" ]; then
    VARIABLE=( job_align_bowtie job_vcf_bowtie )
    METH_FILE_BIG=bowtie2 #bowtie2, aln, mem  
else
    echo "ERROR"
fi

mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/cluster/jobs/

for i in "${!SAMPLES[@]}"
do  
    JOB_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/cluster/jobs/
    FILE_NAME_JOB=${JOB_PATH}"${TISSUE[i]}"_${METHOD}
    echo "${TISSUE[i]}"
    IFS=', ' read -r -a array <<< "${TISSUE[i]}"
    echo "${array[1]}"
    #echo "${hallo[i]}"

    echo '#!/usr/bin/bash' > ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '#SBATCH --job-name=${FILE_NAME_JOB}' >> ${FILE_NAME_JOB}.sh
    echo '#SBATCH --output=${FILE_NAME_JOB}.out' >> ${FILE_NAME_JOB}.sh
    echo '#SBATCH --error=${FILE_NAME_JOB}.err' >> ${FILE_NAME_JOB}.sh
    echo '#SBATCH --time=89:59:59' >> ${FILE_NAME_JOB}.sh
    echo '#SBATCH --cpus-per-task=8' >> ${FILE_NAME_JOB}.sh
    echo '#SBATCH --mem=96gb' >> ${FILE_NAME_JOB}.sh
    echo '#SBATCH --nodes=1' >> ${FILE_NAME_JOB}.sh
    echo '#SBATCH --open-mode=append' >> ${FILE_NAME_JOB}.sh
    echo '#SBATCH --export=NONE' >> ${FILE_NAME_JOB}.sh
    echo '#SBATCH --get-user-env=L' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '# Load packages' >> ${FILE_NAME_JOB}.sh
    echo 'ml BCFtools/1.11-GCCcore-7.3.0' >> ${FILE_NAME_JOB}.sh
    echo 'ml Bowtie2/2.3.4.2-foss-2018b' >> ${FILE_NAME_JOB}.sh
    echo 'ml BWA/0.7.17-GCCcore-7.3.0' >> ${FILE_NAME_JOB}.sh
    echo 'ml cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare' >> ${FILE_NAME_JOB}.sh
    echo 'ml FastQC/0.11.8-Java-11-LTS' >> ${FILE_NAME_JOB}.sh
    echo 'ml GATK/4.1.4.1-Java-8-LTS' >> ${FILE_NAME_JOB}.sh
    echo 'ml Java/11-LTS' >> ${FILE_NAME_JOB}.sh
    echo '#ml libjpeg-turbo/2.0.2-GCCcore-7.3.0 #?' >> ${FILE_NAME_JOB}.sh
    echo 'ml picard/2.20.5-Java-11-LTS' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '#Chrom' >> ${FILE_NAME_JOB}.sh
    echo 'CHROM=${CHROM_CHOOSEN}' >> ${FILE_NAME_JOB}.sh
    echo '# The path to the file of the reference genome that will be used.' >> ${FILE_NAME_JOB}.sh
    echo 'if [ "${CHROM}" != "chrall" ]; then' >> ${FILE_NAME_JOB}.sh
    echo '  PATH_GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/' >> ${FILE_NAME_JOB}.sh
    echo 'else' >> ${FILE_NAME_JOB}.sh
    echo '  PATH_GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/ #chrall' >> ${FILE_NAME_JOB}.sh
    echo 'fi' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo 'GENOOM=${PATH_GENOOM}${CHROM}.fa' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo 'ml SAMtools/1.9-foss-2018b' >> ${FILE_NAME_JOB}.sh
    echo "# When file doesn't exist." >> ${FILE_NAME_JOB}.sh
    echo 'if [ ! -f ${PATH_GENOOM}${CHROM}.dict ]; then' >> ${FILE_NAME_JOB}.sh
    echo '  echo "FILE does not exist."' >> ${FILE_NAME_JOB}.sh
    echo '  # dict file' >> ${FILE_NAME_JOB}.sh
    echo '  gatk CreateSequenceDictionary -R ${GENOOM}' >> ${FILE_NAME_JOB}.sh
    echo '  # index file' >> ${FILE_NAME_JOB}.sh
    echo '  samtools faidx ${GENOOM}' >> ${FILE_NAME_JOB}.sh
    echo 'fi' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '# Array of tissue numbers' >> ${FILE_NAME_JOB}.sh
    echo 'TISSUE_ARR="${array[i]}"' >> ${FILE_NAME_JOB}.sh
    echo '# Array of sample numbers' >> ${FILE_NAME_JOB}.sh
    echo 'SAMPLE="${SAMPLES[i]}"' >> ${FILE_NAME_JOB}.sh #TODO TODO TODO
    echo 'array3=${SAMPLES}' >> ${FILE_NAME_JOB}.sh #TODO TODO TODO
    echo '' >> ${FILE_NAME_JOB}.sh
    echo 'METHOD=${METHOD_BIG} #bowtie,   bwa_aln, bwa_mem' >> ${FILE_NAME_JOB}.sh
    echo 'METH_FILE=${METH_FILE_BIG} #bowtie2, aln, mem' >> ${FILE_NAME_JOB}.sh
    echo 'TYPE_ALN=mutect_${METHOD}' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo 'GENERAL_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/' >> ${FILE_NAME_JOB}.sh
    echo 'SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/cluster/' >> ${FILE_NAME_JOB}.sh
    echo '# Path to the file of the panel of normals (PoN)' >> ${FILE_NAME_JOB}.sh
    echo 'PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/merge_somatic-b37_Mutect2-WGS-panel-b37.vcf' >> ${FILE_NAME_JOB}.sh
    echo '# Path to the file of the germline resource (GR)' >> ${FILE_NAME_JOB}.sh
    echo 'GR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GermlineResource/merge_af-only-gnomad.raw.sites.vcf' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '# for python file parameters: automatic_script_ob.py' >> ${FILE_NAME_JOB}.sh
    echo '# Number of tumors you want to combine while running Mutect2.' >> ${FILE_NAME_JOB}.sh
    echo 'NUMBER_OF_TUMORS_py=1' >> ${FILE_NAME_JOB}.sh
    echo '# Number of hc you want to combine while running Mutect2' >> ${FILE_NAME_JOB}.sh
    echo 'NUMBER_OF_HC_py=1' >> ${FILE_NAME_JOB}.sh
    echo '# What type of tumor you want to have ("both", "tFL" or "FL")' >> ${FILE_NAME_JOB}.sh
    echo 'TYPE_SAMPLE_py='both'' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '# load conda and activate to execute python script' >> ${FILE_NAME_JOB}.sh
    echo 'ml Anaconda3/5.3.0' >> ${FILE_NAME_JOB}.sh
    echo 'source activate stage' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '# RUN: file_prep.sh' >> ${FILE_NAME_JOB}.sh
    echo '#source ${SCRIPT_PATH}file_prep.sh' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo 'echo ${METHOD}' >> ${FILE_NAME_JOB}.sh
    echo '# RUN: "${VARIABLE[0]}".sh' >> ${FILE_NAME_JOB}.sh
    echo 'echo "${VARIABLE[0]}"' >> ${FILE_NAME_JOB}.sh
    echo 'source ${SCRIPT_PATH}"${VARIABLE[0]}".sh' >> ${FILE_NAME_JOB}.sh
    echo 'for i in "${array3[@]}"' >> ${FILE_NAME_JOB}.sh ## TODO TODO
    echo 'do' >> ${FILE_NAME_JOB}.sh
    echo '  mkdir -p ${GENERAL_PATH}"${i}"/${CHROM}/mutect_${METHOD}/' >> ${FILE_NAME_JOB}.sh
    echo 'done' >> ${FILE_NAME_JOB}.sh
    echo '# RUN: automatic_script_ob.py' >> ${FILE_NAME_JOB}.sh
    echo "echo 'automatic_script_ob'" >> ${FILE_NAME_JOB}.sh
    echo 'python3 ${SCRIPT_PATH}automatic_script_ob.py ${GENERAL_PATH} ${NUMBER_OF_TUMORS_py} ${NUMBER_OF_HC_py} ${TYPE_SAMPLE_py} ${METHOD} ${METH_FILE} ${CHROM}' >> ${FILE_NAME_JOB}.sh
    echo '# RUN: change_sample_name.sh' >> ${FILE_NAME_JOB}.sh
    echo "echo 'change_sample_name'" >> ${FILE_NAME_JOB}.sh
    echo 'source ${SCRIPT_PATH}change_sample_name.sh' >> ${FILE_NAME_JOB}.sh
    echo '# RUN: "${VARIABLE[0]}".sh' >> ${FILE_NAME_JOB}.sh
    echo '# mutect2' >> ${FILE_NAME_JOB}.sh
    echo 'echo "${VARIABLE[0]}"' >> ${FILE_NAME_JOB}.sh
    echo 'source ${SCRIPT_PATH}"${VARIABLE[0]}".sh' >> ${FILE_NAME_JOB}.sh
    echo '# RUN: vcf_compare_auto.sh' >> ${FILE_NAME_JOB}.sh
    echo "echo 'vcf_compare_auto'" >> ${FILE_NAME_JOB}.sh
    echo 'source ${SCRIPT_PATH}vcf_compare_auto.sh' >> ${FILE_NAME_JOB}.sh
    echo '# RUN: vcf_merge_auto.sh' >> ${FILE_NAME_JOB}.sh
    echo '# Merge files' >> ${FILE_NAME_JOB}.sh
    echo 'COMP_TYPE=mutect2 # manual, mutect2' >> ${FILE_NAME_JOB}.sh
    echo "echo 'vcf_merge_auto1'" >> ${FILE_NAME_JOB}.sh
    echo 'source ${SCRIPT_PATH}vcf_merge_auto.sh' >> ${FILE_NAME_JOB}.sh
    echo 'COMP_TYPE=manual # manual, mutect2' >> ${FILE_NAME_JOB}.sh
    echo "echo 'vcf_merge_auto2'" >> ${FILE_NAME_JOB}.sh
    echo 'source ${SCRIPT_PATH}vcf_merge_auto.sh' >> ${FILE_NAME_JOB}.sh
    echo '# RUN: annotate.sh' >> ${FILE_NAME_JOB}.sh
    echo '# filteren van dbSNP (annoteren)' >> ${FILE_NAME_JOB}.sh
    echo "echo 'annotate'" >> ${FILE_NAME_JOB}.sh
    echo 'source ${SCRIPT_PATH}annotate.sh' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo '' >> ${FILE_NAME_JOB}.sh
    echo "echo 'THE END END END END'" >> ${FILE_NAME_JOB}.sh











done