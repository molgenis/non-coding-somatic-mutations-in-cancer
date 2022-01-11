#!/usr/bin/bash

NUM=5042
FILE_1=SS600${NUM}_chr22.1.fastq
FILE_2=SS600${NUM}_chr22.2.fastq
OUTPUT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NUM}/QC/

cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NUM}/

ml FastQC/0.11.8-Java-11-LTS
ml cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare

echo "BEGIN"
# fastqc
# https://wiki.bits.vib.be/index.php/Linux_command_line#Automating_FASTQC_analyses
fastqc -f fastq -o ${OUTPUT} ${FILE_1} ${FILE_2}
# fastqc -f bam -o ${OUTPUT} SS6005044_chr22.sort.bam


# trimming pair-end
# https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads
cutadapt -o QC/cutadapt/TRIM_${FILE_1} -p QC/cutadapt/TRIM_${FILE_2} ${FILE_1} ${FILE_2} -q 20 --minimum-length 30 --pair-filter=any
fastqc -f fastq -o ${OUTPUT}/cutadapt/ QC/cutadapt/TRIM_${FILE_1} QC/cutadapt/TRIM_${FILE_2}

# BBDuk
#https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
#Length filtering:
#This will discard reads shorter than 100bp after trimming to Q10. 
#Alternatively, using “mlf=50” (minlengthfraction=50) would discard reads under 50% of their original length after trimming. 
#Either of these flags apply to any trim operation, whether force-trim (ftr, ftl, ftm), quality-trimming (qtrim), or kmer-trimming (ktrim). 
#“mlf” compares the final length after all trim operations to the initial length before any trim operations.

#Quality trimming:
# This will quality-trim to Q10 using the Phred algorithm, which is more accurate than naive trimming. 
# “qtrim=r” means it will trim the right side only; you can alternatively set “qtrim=l” for left or “qtrim=rl” for both. 
# If quality trimming is enabled, it happens after all kmer-based operations.
bbduk.sh -Xmx1g in1=${FILE_1} in2=${FILE_2} out1=QC/bbduk/clean_20_${FILE_1} out2=QC/bbduk/clean_20_${FILE_2} qtrim=rl trimq=20 minlen=30
fastqc -f fastq -o ${OUTPUT}/bbduk/ QC/bbduk/clean_20_${FILE_1} QC/bbduk/clean_20_${FILE_2}
# Quality filtering:
# This will discard reads with average quality below 20. If quality-trimming is enabled, the average quality will be calculated on the trimmed read.
# NIET TRIMMEN
bbduk.sh -Xmx1g in1=${FILE_1} in2=${FILE_2} out1=QC/bbduk/clean_qf20_${FILE_1} out2=QC/bbduk/clean_qf20__${FILE_2} maq=20
fastqc -f fastq -o ${OUTPUT}/bbduk/ QC/bbduk/clean_qf20_${FILE_1} QC/bbduk/clean_qf20__${FILE_2}


#http://hannonlab.cshl.edu/fastx_toolkit/
#conda install -c bioconda fastx_toolkit
# FASTQ Quality Filter: Filters sequences based on quality
#fastq_quality_filter


echo "EIND"
