#!/usr/bin/bash

#SBATCH --job-name=pipelineS3
#SBATCH --output=pipelineS3.out
#SBATCH --error=pipelineS3.err
#SBATCH --time=89:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

# Load packages
ml BCFtools/1.11-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml BWA/0.7.17-GCCcore-7.3.0
ml cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare
ml FastQC/0.11.8-Java-11-LTS
ml GATK/4.1.4.1-Java-8-LTS
ml Java/11-LTS
####ml libjpeg-turbo/2.0.2-GCCcore-7.3.0 #?
ml picard/2.20.5-Java-11-LTS


# Chosen chromosome
CHROM=chr21
# Array of tissue numbers
array=(  4094  4099  4104  4109  4113  4114  4118  4119  4123  4124  4129  4128  4133  4134  4139  4138  5041  5043  5042  5044  ) 
# Array of sample numbers
array2=(  S1  S1  S1  S2  S2  S2  S2  S2  S3  S3  S3  S4  S4  S4  S5  S5  S5  S6  S6  S6 ) 
# FOR job_vcf_aln.sh............
array3=( S1  S2  S3  S4  S5  S6 ) 
# Chosen alignment method
METHOD=bwa_aln #bowtie, bwa_aln, bwa_mem
METH_FILE=aln  #bowtie2, aln, mem
TYPE_ALN=mutect_${METHOD}
# The global path to where the data will be and is placed
GENERAL_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/
# Path to the scripts called in this script and other scripts
SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/cluster/
# Path to the file of the panel of normals (PoN)
PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/merge_somatic-b37_Mutect2-WGS-panel-b37.vcf
# Path to the file of the germline resource (GR)
GR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GermlineResource/merge_af-only-gnomad.raw.sites.vcf

# For python file parameters: automatic_script_ob.py
# Number of tumors you want to combine while running Mutect2.
NUMBER_OF_TUMORS_py=1
# Number of hc you want to combine while running Mutect2
NUMBER_OF_HC_py=1
# What type of tumor you want to have ("both", "tFL" or "FL")
TYPE_SAMPLE_py='both'

# The path to the file of the reference genome that will be used.
if [ "${CHROM}" != "chrall" ]; then
    PATH_GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/
else
    PATH_GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/ #chrall
fi

# Path to the genome file
GENOOM=${PATH_GENOOM}${CHROM}.fa

# Load samtools
ml SAMtools/1.9-foss-2018b
# When file doesn't exist.
if [ ! -f ${PATH_GENOOM}${CHROM}.dict ]; then
    echo "FILE does not exist."
    # dict file
    gatk CreateSequenceDictionary -R ${GENOOM}
    # index file
    samtools faidx ${GENOOM}
fi

# load conda and activate to execute python script
ml Anaconda3/5.3.0
source activate stage

# RUN: file_prep.sh
#source ${SCRIPT_PATH}file_prep.sh

# Lexicographic (greater than, less than) comparison.
if [ "${METHOD}" == "bwa_aln" ]; then
    echo ${METHOD}
    # RUN: job_align_aln.sh
    source ${SCRIPT_PATH}job_align_aln.sh
    # Make different folders
    for i in "${array3[@]}"
    do 
        mkdir -p ${GENERAL_PATH}"${i}"/${CHROM}/mutect_${METHOD}/
    done
    # RUN: automatic_script_ob.py
    python3 ${SCRIPT_PATH}automatic_script_ob.py ${GENERAL_PATH} ${NUMBER_OF_TUMORS_py} ${NUMBER_OF_HC_py} ${TYPE_SAMPLE_py} ${METHOD} ${METH_FILE} ${CHROM}
    # RUN: change_sample_name.sh
    source ${SCRIPT_PATH}change_sample_name.sh
    # RUN: job_vcf_aln.sh
    # DUURT ERG LANG!
    source ${SCRIPT_PATH}job_vcf_aln.sh    
    # RUN: vcf_compare_auto.sh
    source ${SCRIPT_PATH}vcf_compare_auto.sh
    
    
    #COMP_TYPE=mutect2 # manual, mutect2
    # source ${SCRIPT_PATH}vcf_merge_auto.sh
    # source ${SCRIPT_PATH}annotate.sh
    # source ${SCRIPT_PATH}FORMAT/write_df.sh
    # source ${SCRIPT_PATH}FORMAT/make_plots_loop.sh
    # ml R/4.0.3-foss-2018b-bare
    # make_file_for_chromosome_plots.py
    # split.py
    # karyoploteR_plots.R
    # chromplot_plots.R
    

elif [ "${METHOD}" == "bwa_mem" ]; then
    echo ${METHOD}
    
elif [ "${METHOD}" == "bowtie" ]; then
    echo ${METHOD}
    
else
    echo "ERROR"
fi

echo 'THE END END END END'