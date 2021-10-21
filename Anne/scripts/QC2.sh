#!/usr/bin/bash

NUM=${1}
FILE_1=SS600${NUM}_chr22.1.fastq
FILE_2=SS600${NUM}_chr22.2.fastq
OUTPUT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NUM}/QC/
LEN=${2}

echo "${LEN}"

mkdir -p /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NUM}/QC/cutadapt/

cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NUM}/

ml FastQC/0.11.8-Java-11-LTS
ml cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare

echo "BEGIN"
# fastqc
# https://wiki.bits.vib.be/index.php/Linux_command_line#Automating_FASTQC_analyses
#fastqc -f fastq -o ${OUTPUT} ${FILE_1} ${FILE_2}
# fastqc -f bam -o ${OUTPUT} SS6005044_chr22.sort.bam


# trimming pair-end
# https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads
cutadapt -o QC/cutadapt/TRIM_${LEN}_${FILE_1} -p QC/cutadapt/TRIM_${LEN}_${FILE_2} ${FILE_1} ${FILE_2} -q 20 --minimum-length ${LEN} 1> QC/cutadapt/summary_for_TRIM_${LEN}_${NUM}.txt #--pair-filter=any
fastqc -f fastq -o ${OUTPUT}/cutadapt/ QC/cutadapt/TRIM_${LEN}_${FILE_1} QC/cutadapt/TRIM_${LEN}_${FILE_2}

#http://hannonlab.cshl.edu/fastx_toolkit/
#http://hannonlab.cshl.edu/fastx_toolkit/install_ubuntu.txt
#http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_quality_filter_usage
#conda install -c bioconda fastx_toolkit
# FASTQ Quality Filter: Filters sequences based on quality
#fastq_quality_filter

#fastqc -f fastq -o ${OUTPUT} ${FILE_1} ${FILE_2}




echo "EIND"
