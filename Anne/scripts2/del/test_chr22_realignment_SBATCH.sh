#!/usr/bin/env bash

#SBATCH --job-name=testchr22alignment
#SBATCH --output=testchr22alignment.out
#SBATCH --error=testchr22alignment.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

#ml SAMtools
#ml HTSlib
#ml picard
cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292

ml SAMtools/1.9-foss-2018b
ml GATK/4.1.4.1-Java-8-LTS
ml BWA/0.7.17-GCCcore-7.3.0
ml picard/2.20.5-Java-11-LTS

#whateverotherlibs
FILE="SS6005043"

samtools view -h ${FILE}.sorted.bam chr22 > ${FILE}.chr22.sam
samtools view -bS ${FILE}.chr22.sam > ${FILE}.chr22.bam
samtools fastq ${FILE}.chr22.bam > ${FILE}.chr22.fastq
gzip ${FILE}.chr22.fastq
bwa index -a is /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa
java -jar ${EBROOTPICARD}/picard.jar CreateSequenceDictionary R= /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa O= /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.dict
samtools faidx /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa
bwa aln /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa ${FILE}.chr22.fastq.gz -n 0.1 -l 1000 > ${FILE}.chr22.sai
bwa samse /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa ${FILE}.chr22.sai ${FILE}.chr22.fastq.gz -f ${FILE}.chr22.sam
samtools view -Sb ${FILE}.chr22.sam > ${FILE}.chr22.bam
samtools view ${FILE}.chr22.bam | less -S
samtools view -H ${FILE}.chr22.bam
samtools sort ${FILE}.chr22.bam -o ${FILE}.chr22.sort.bam
samtools index ${FILE}.chr22.sort.bam
java -jar ${EBROOTPICARD}/picard.jar  AddOrReplaceReadGroups INPUT= ${FILE}.chr22.sort.bam OUTPUT= ${FILE}.chr22.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT
samtools index ${FILE}.chr22.RG.bam
java -jar ${EBROOTPICARD}/picard.jar  MarkDuplicates I= ${FILE}.chr22.RG.bam O= ${FILE}.chr22.DR.bam M=output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> logFile.log
samtools index ${FILE}.chr22.DR.bam

        
