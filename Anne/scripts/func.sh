#!/usr/bin/bash

#FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/SS6004099.sorted.bam
#TYPE=tumor #tumor or normal

TUMOR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/
NORMAL=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/

PATH_REF=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/
REF=chr22.fa

CHR_DICT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.dict
CHR=chr22

VCF_PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/
VCF_NAME=somatic-b37_Mutect2-WGS-panel-b37.vcf
PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/merge_somatic-b37_Mutect2-WGS-panel-b37.vcf

ml Anaconda3/5.3.0
ml SAMtools/1.9-foss-2018b
ml GATK/4.1.4.1-Java-8-LTS
ml BCFtools/1.11-GCCcore-7.3.0
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS




change_sample_name() {
    samtools view -H ${1}.${2}.bam  | sed "s/SM:[^\t]*/SM:${1}/g" | samtools reheader - ${1}.${2}.bam > SN_${1}.${2}.bam
    samtools index SN_${1}.${2}.bam
}

Mutect2_one_sample(){
    # standard
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${3}.bam -O unfiltered_${1}_${CHR}_${3}_somatic.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${CHR}_${3}_somatic.vcf.gz -O filtered_${1}_${CHR}_${3}_somatic.vcf.gz
    # sample name
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${3}.bam -${2} ${1}.${3}.bam -O unfiltered_${1}_${CHR}_${3}_somatic_name.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${CHR}_${3}_somatic_name.vcf.gz -O filtered_${1}_${CHR}_${3}_somatic_name.vcf.gz
    # sample name and PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${3}.bam -${2} ${1}.${3}.bam -O unfiltered_${1}_${CHR}_${3}_somatic_name_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${CHR}_${3}_somatic_name_PON.vcf.gz -O filtered_${1}_${CHR}_${3}_somatic_name_PON.vcf.gz
    # PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${3}.bam -O unfiltered_${1}_${CHR}_${3}_somatic_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${CHR}_${3}_somatic_PON.vcf.gz -O filtered_${1}_${CHR}_${3}_somatic_PON.vcf.gz   
}




#Filter_with_hand() {
#}

cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/

printf "####change_sample_name\n"
change_sample_name aln-pe RG #1=file
change_sample_name aln-pe DR #1=file

printf "\n####Mutect2_one_sample\n"
Mutect2_one_sample aln-pe tumor RG #1=file 2=type
Mutect2_one_sample aln-pe tumor DR #1=file 2=type





