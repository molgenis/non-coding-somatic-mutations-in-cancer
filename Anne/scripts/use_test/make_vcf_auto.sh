#!/usr/bin/bash

PATH_REF=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/
REF=chr22.fa

#CHR_DICT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.dict
CHR=chr22

#VCF_PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/
#VCF_NAME=somatic-b37_Mutect2-WGS-panel-b37.vcf
PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/merge_somatic-b37_Mutect2-WGS-panel-b37.vcf
GR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GermlineResource/merge_af-only-gnomad.raw.sites.vcf

ml Anaconda3/5.3.0
ml SAMtools/1.9-foss-2018b
ml GATK/4.1.4.1-Java-8-LTS
ml BCFtools/1.11-GCCcore-7.3.0
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS

# Function that changed the sample names
change_sample_name() {
    samtools view -H ${2}${1}.bam  | sed "s/SM:[^\t]*/SM:${1}/g" | samtools reheader - ${2}${1}.bam > ${2}SN_${1}.bam
    samtools index ${2}SN_${1}.bam
}

echo_test(){
    echo ${1}
    echo ${2}
    echo ${3}
}

# bam --> VCF
mutect2_vcf() {
    mkdir -p ${1}
    # sample name
    echo "####sample name\n"
    gatk Mutect2 -R ${PATH_REF}${REF} ${2} -O ${3}_somatic_unfiltered.vcf.gz --native-pair-hmm-threads 8
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}_somatic_unfiltered.vcf.gz -O ${3}_somatic_filtered.vcf.gz
    # PON
    echo "####PON\n"
    gatk Mutect2 -R ${PATH_REF}${REF} ${2} -O ${3}_somatic_unfiltered_PON.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}_somatic_unfiltered_PON.vcf.gz -O ${3}_somatic_filtered_PON.vcf.gz
    #GR
    echo "####GR\n"
    gatk Mutect2 -R ${PATH_REF}${REF} ${2} -O ${3}_somatic_unfiltered_GERM.vcf.gz --native-pair-hmm-threads 8 --germline-resource ${GR}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} --V ${3}_somatic_unfiltered_GERM.vcf.gz -O ${3}_somatic_filtered_GERM.vcf.gz
    #PON and GR
    echo "####PON and GR\n"
    gatk Mutect2 -R ${PATH_REF}${REF} ${2} -O ${3}_somatic_unfiltered_PON_GERM.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8 --germline-resource ${GR}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} --V ${3}_somatic_unfiltered_PON_GERM.vcf.gz -O ${3}_somatic_filtered_PON_GERM.vcf.gz
}

