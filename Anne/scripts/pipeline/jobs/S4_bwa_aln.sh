#!/usr/bin/bash

#SBATCH --job-name=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/pipeline/jobs/S4_bwa_aln
#SBATCH --output=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/pipeline/jobs/S4_bwa_aln.out
#SBATCH --error=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/pipeline/jobs/S4_bwa_aln.err
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
#ml libjpeg-turbo/2.0.2-GCCcore-7.3.0 #?
ml picard/2.20.5-Java-11-LTS

#Chrom
CHROM=chr21
# The path to the file of the reference genome that will be used.
if [ "${CHROM}" != "chrall" ]; then
  PATH_GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/
else
  PATH_GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/ #chrall
fi

GENOOM=${PATH_GENOOM}${CHROM}.fa

ml SAMtools/1.9-foss-2018b
# When file doesn't exist.
if [ ! -f ${PATH_GENOOM}${CHROM}.dict ]; then
  echo 'FILE does not exist.'
  # dict file
  gatk CreateSequenceDictionary -R ${GENOOM}
  # index file
  samtools faidx ${GENOOM}
fi

# Array of tissue numbers
TISSUE_ARR=(4128 4133 4134)
# Sample numbers
SAMPLE=S4

METHOD=bwa_aln #bowtie,   bwa_aln, bwa_mem
METH_FILE=aln #bowtie2, aln, mem
TYPE_ALN=mutect_${METHOD}

GENERAL_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/
SCRIPT_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/pipeline/
# Path to the file of the panel of normals (PoN)
PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/merge_somatic-b37_Mutect2-WGS-panel-b37.vcf
# Path to the file of the germline resource (GR)
GR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GermlineResource/merge_af-only-gnomad.raw.sites.vcf

# for python file parameters: automatic_script_ob.py
# Number of tumors you want to combine while running Mutect2.
NUMBER_OF_TUMORS_py=1
# Number of hc you want to combine while running Mutect2
NUMBER_OF_HC_py=1
# What type of tumor you want to have ('both', 'tFL' or 'FL')
TYPE_SAMPLE_py='both'

# load conda and activate to execute python script
ml Anaconda3/5.3.0
source activate stage

# RUN: file_prep.sh
#source ${SCRIPT_PATH}file_prep.sh

echo ${METHOD}
# RUN: align.sh
echo 'align.sh'
source ${SCRIPT_PATH}align.sh
mkdir -p ${GENERAL_PATH}${SAMPLE}/${CHROM}/mutect_${METHOD}/
# RUN: make_file_vc.py
echo 'make_file_vc.py'
python3 ${SCRIPT_PATH}make_file_vc.py ${GENERAL_PATH} ${NUMBER_OF_TUMORS_py} ${NUMBER_OF_HC_py} ${TYPE_SAMPLE_py} ${METHOD} ${METH_FILE} ${CHROM}
# RUN: change_sample_name.sh
echo 'change_sample_name'
source ${SCRIPT_PATH}change_sample_name.sh
# RUN: variant_calling_arguments.sh
# mutect2
echo 'variant_calling_arguments.sh'
source ${SCRIPT_PATH}variant_calling_arguments.sh
# RUN: compare_vcf.sh
echo 'compare_vcf.sh'
source ${SCRIPT_PATH}compare_vcf.sh
# RUN: merge_vcf.sh
# Merge files
COMP_TYPE=mutect2 # manual, mutect2
echo 'merge_vcf1'
source ${SCRIPT_PATH}merge_vcf.sh
COMP_TYPE=manual # manual, mutect2
echo 'vmerge_vcf2'
source ${SCRIPT_PATH}merge_vcf.sh
# RUN: dbSNP_annotate.sh
# filteren van dbSNP (annoteren)
echo 'dbSNP_annotate'
source ${SCRIPT_PATH}dbSNP_annotate.sh



echo 'THE END END END END'
