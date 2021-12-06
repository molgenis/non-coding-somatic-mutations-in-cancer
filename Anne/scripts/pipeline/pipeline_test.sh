#!/usr/bin/bash

# Load packages
ml Anaconda3/5.3.0
ml BCFtools/1.11-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml BWA/0.7.17-GCCcore-7.3.0
ml cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare
ml FastQC/0.11.8-Java-11-LTS
ml GATK/4.1.4.1-Java-8-LTS
ml picard/2.20.5-Java-11-LTS
ml SAMtools/1.9-foss-2018b



CHROM=chr21
# The path to the file of the reference genome that will be used.
if [ "${CHROM}" != "chrall" ]; then
    GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHROM}.fa
else
    GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/${CHROM}.fa #chrall
fi

# dict file
gatk CreateSequenceDictionary -R ${GENOOM}
# index file
samtools faidx ${GENOOM}


# Array of tissue numbers
array=( 4094 4099) #  4101  ) #4109  4113  4114  4118  4119  4123  4124  4129  4128  4133  4134  4139  4138  5041  5043  5042  5044 )
# Array of sample numbers
array2=( S1  S1 )#  S1 ) # S2  S2  S2  S2  S2  S3  S3  S3  S4  S4  S4  S5  S5  S5  S6  S6  S6 )

METHOD=bwa_aln #bowtie,   bwa_aln, bwa_mem
METH_FILE=aln #bowtie2, aln, mem

GENERAL_PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/

# for python file parameters: automatic_script_ob.py
# Number of tumors you want to combine while running Mutect2.
NUMBER_OF_TUMORS_py=1
# Number of hc you want to combine while running Mutect2
NUMBER_OF_HC_py=1
# What type of tumor you want to have ("both", "tFL" or "FL")
TYPE_SAMPLE_py='both'




# RUN: file_prep.sh
# source ${GENERAL_PATH}file_prep.sh

# Lexicographic (greater than, less than) comparison.
if [ "${METHOD}" == "bwa_aln" ]; then
    echo ${METHOD}
    # RUN: job_align_aln.sh
    # source ${GENERAL_PATH}job_align_aln.sh
    # RUN: automatic_script_ob.py
    # python3 ${GENERAL_PATH}automatic_script_ob.py ${GENERAL_PATH} ${NUMBER_OF_TUMORS_py} ${NUMBER_OF_HC_py} ${TYPE_SAMPLE_py} ${METHOD} ${METH_FILE} ${CHROM}
    # RUN: job_vcf_aln.sh
    # source ${GENERAL_PATH}job_vcf_aln.sh

elif [ "${METHOD}" == "bwa_mem" ]; then
    echo ${METHOD}
    
elif [ "${METHOD}" == "bowtie" ]; then
    echo ${METHOD}
    
else
    echo "ERROR"
fi

echo 'THE END END END END'