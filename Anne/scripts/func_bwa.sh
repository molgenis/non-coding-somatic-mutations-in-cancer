#!/usr/bin/bash

#FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/SS6004099.sorted.bam
#TYPE=tumor #tumor or normal

PATH_REF=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/
REF=chr22.fa

CHR_DICT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.dict
CHR=chr22

VCF_PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/
VCF_NAME=somatic-b37_Mutect2-WGS-panel-b37.vcf
PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/merge_somatic-b37_Mutect2-WGS-panel-b37.vcf

dir=5042_5044

ml Anaconda3/5.3.0
ml SAMtools/1.9-foss-2018b
ml GATK/4.1.4.1-Java-8-LTS
ml BCFtools/1.11-GCCcore-7.3.0
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS




change_sample_name() {
    samtools view -H ${2}${1}.bam  | sed "s/SM:[^\t]*/SM:${1}/g" | samtools reheader - ${2}${1}.bam > ${2}SN_${1}.bam
    samtools index ${2}SN_${1}.bam
}

Mutect2_one_sample(){
    # standard
    printf "####standard\n"
    gatk Mutect2 -R ${PATH_REF}${REF} -I ${3}SN_${1}.bam -O ${3}unfiltered_${1}_${CHR}_somatic.vcf.gz --native-pair-hmm-threads 8
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${1}_${CHR}_somatic.vcf.gz -O ${3}filtered_${1}_${CHR}_somatic.vcf.gz
    # sample name
    printf "####sample name\n"
    gatk Mutect2 -R ${PATH_REF}${REF} -I ${3}SN_${1}.bam -${2} ${1}.bam -O ${3}unfiltered_${1}_${CHR}_somatic_name.vcf.gz --native-pair-hmm-threads 8
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${1}_${CHR}_somatic_name.vcf.gz -O ${3}filtered_${1}_${CHR}_somatic_name.vcf.gz
    # sample name and PON
    printf "####sample name and PON\n"
    gatk Mutect2 -R ${PATH_REF}${REF} -I ${3}SN_${1}.bam -${2} ${1}.bam -O ${3}unfiltered_${1}_${CHR}_somatic_name_PON.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${1}_${CHR}_somatic_name_PON.vcf.gz -O ${3}filtered_${1}_${CHR}_somatic_name_PON.vcf.gz
    # PON
    printf "####PON\n"
    gatk Mutect2 -R ${PATH_REF}${REF} -I ${3}SN_${1}.bam -O ${3}unfiltered_${1}_${CHR}_somatic_PON.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${1}_${CHR}_somatic_PON.vcf.gz -O ${3}filtered_${1}_${CHR}_somatic_PON.vcf.gz   
}

Mutect2_two_samples() {
    # standard
    gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -O ${OUTPUT_DICT}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${OUTPUT_DICT}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic.vcf.gz -O ${OUTPUT_DICT}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic.vcf.gz
    # sample name
    gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -tumor ${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -normal ${2}.bam -O ${OUTPUT_DICT}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz 
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${OUTPUT_DICT}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz -O ${OUTPUT_DICT}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz
    # sample name and PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -tumor ${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -normal ${2}.bam -O ${OUTPUT_DICT}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${OUTPUT_DICT}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz -O ${OUTPUT_DICT}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz
    # PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -O ${OUTPUT_DICT}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${OUTPUT_DICT}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz -O ${OUTPUT_DICT}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz
}



cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${dir}

TUMOR=5044
NORMAL=5042

PATH_DICT_one=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${TUMOR}/bwa_aln/
PATH_DICT_two=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NORMAL}/bwa_aln/
OUTPUT_DICT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/${NORMAL}_${TUMOR}/bwa_aln/



printf "####change_sample_name\n"
change_sample_name SS600${TUMOR}_${CHR}_mapped_1_2.bwa.RG ${PATH_DICT_one} #1=file
change_sample_name SS600${TUMOR}_${CHR}_mapped_1_2.bwa.DR ${PATH_DICT_one} #1=file
change_sample_name SS600${NORMAL}_${CHR}_mapped_1_2.bwa.RG ${PATH_DICT_two} #1=file
change_sample_name SS600${NORMAL}_${CHR}_mapped_1_2.bwa.DR ${PATH_DICT_two} #1=file


printf "\n####Mutect2_one_sample\n"
#Mutect2_one_sample SS600${TUMOR}_${CHR}_mapped_1_2.bwa.RG tumor ${PATH_DICT_one} #1=file 2=type
#Mutect2_one_sample SS600${TUMOR}_${CHR}_mapped_1_2.bwa.DR tumor ${PATH_DICT_one} #1=file 2=type
Mutect2_one_sample SS600${NORMAL}_${CHR}_mapped_1_2.bwa.RG normal ${PATH_DICT_two} #1=file 2=type
Mutect2_one_sample SS600${NORMAL}_${CHR}_mapped_1_2.bwa.DR normal ${PATH_DICT_two} #1=file 2=type

printf "\n####Mutect2_two_sample\n"
printf "###########5044-5042\n"
Mutect2_two_samples SS600${TUMOR}_${CHR}_mapped_1_2.bwa.RG SS600${NORMAL}_${CHR}_mapped_1_2.bwa.RG #1=tumor 2=normal
Mutect2_two_samples SS600${TUMOR}_${CHR}_mapped_1_2.bwa.DR SS600${NORMAL}_${CHR}_mapped_1_2.bwa.DR  #1=tumor 2=normal




