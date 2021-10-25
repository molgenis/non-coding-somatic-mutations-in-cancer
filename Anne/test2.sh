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

# # Function that changed the sample names
# change_sample_name() {
#     samtools view -H ${2}${1}.bam  | sed "s/SM:[^\t]*/SM:${1}/g" | samtools reheader - ${2}${1}.bam > ${2}SN_${1}.bam
#     samtools index ${2}SN_${1}.bam
# }

# # bam --> VCF, for one file
# Mutect2_one_sample(){
#     # standard
#     #echo "####standard\n"
#     #gatk Mutect2 -R ${PATH_REF}${REF} -I ${3}SN_${1}.bam -O ${3}unfiltered_${1}_${CHR}_somatic.vcf.gz --native-pair-hmm-threads 8
#     #gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${1}_${CHR}_somatic.vcf.gz -O ${3}filtered_${1}_${CHR}_somatic.vcf.gz
#     # sample name
#     #echo "####sample name\n"
#     #gatk Mutect2 -R ${PATH_REF}${REF} -I ${3}SN_${1}.bam -${2} ${1}.bam -O ${3}unfiltered_${1}_${CHR}_somatic_name.vcf.gz --native-pair-hmm-threads 8
#     #gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${1}_${CHR}_somatic_name.vcf.gz -O ${3}filtered_${1}_${CHR}_somatic_name.vcf.gz
#     # sample name and PON
#     #echo "####sample name and PON\n"
#     #gatk Mutect2 -R ${PATH_REF}${REF} -I ${3}SN_${1}.bam -${2} ${1}.bam -O ${3}unfiltered_${1}_${CHR}_somatic_name_PON.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8
#     #gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${1}_${CHR}_somatic_name_PON.vcf.gz -O ${3}filtered_${1}_${CHR}_somatic_name_PON.vcf.gz
#     # PON
#     #echo "####PON\n"
#     #gatk Mutect2 -R ${PATH_REF}${REF} -I ${3}SN_${1}.bam -O ${3}unfiltered_${1}_${CHR}_somatic_PON.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8
#     #gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${1}_${CHR}_somatic_PON.vcf.gz -O ${3}filtered_${1}_${CHR}_somatic_PON.vcf.gz
#      sample name and PON and GR
#     echo "####sample name and PON and GR\n"
#     gatk Mutect2 -R ${PATH_REF}${REF} -I ${3}SN_${1}.bam -${2} ${1}.bam -O ${3}unfiltered_${1}_${CHR}_somatic_name_PON_GR.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8 --germline-resource ${GR}
#     gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${1}_${CHR}_somatic_name_PON_GR.vcf.gz -O ${3}filtered_${1}_${CHR}_somatic_name_PON_GR.vcf.gz
# }

# # bam --> VCF, for two files (one tumor file and one normal file)
# Mutect2_two_samples() {
#     # standard
#     echo "####standard\n"
#     gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -O ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic.vcf.gz --native-pair-hmm-threads 8
#     gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic.vcf.gz -O ${3}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic.vcf.gz
#     # sample name
#     echo "####sample name\n"
#     gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -tumor ${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -normal ${2}.bam -O ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz --native-pair-hmm-threads 8
#     gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz -O ${3}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz
#     # sample name and PON
#     echo "####sample name and PON\n"
#     gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -tumor ${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -normal ${2}.bam -O ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8
#     gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz -O ${3}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz
#     # PON
#     echo "####PON\n"
#     gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -O ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8
#     gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz -O ${3}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz
#     # sample name and PON and GR
#     echo "####sample name and PON and GR\n"
#     gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -tumor ${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -normal ${2}.bam -O ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON_GR.vcf.gz --panel-of-normals ${PON} --germline-resource ${GR} --native-pair-hmm-threads 8
#     gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON_GR.vcf.gz -O ${3}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON_GR.vcf.gz
# }


# bam --> VCF, for two files (one tumor file and one normal file)
test1() {
    # standard
    echo "####standard\n"
    ${1}
    gatk Mutect2 -R ${PATH_REF}${REF} ${2} --native-pair-hmm-threads 8
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}_somatic.vcf.gz -O ${4}_somatic.vcf.gz
    # # sample name
    # echo "####sample name\n"
    # gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -tumor ${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -normal ${2}.bam -O ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz --native-pair-hmm-threads 8
    # gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz -O ${3}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz
    # # sample name and PON
    # echo "####sample name and PON\n"
    # gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -tumor ${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -normal ${2}.bam -O ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8
    # gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz -O ${3}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz
    # # PON
    # echo "####PON\n"
    # gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -O ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8
    # gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz -O ${3}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz
    # # sample name and PON and GR
    # echo "####sample name and PON and GR\n"
    # gatk Mutect2 -R ${PATH_REF}${REF} -I ${PATH_DICT_one}SN_${1}.bam -tumor ${1}.bam -I ${PATH_DICT_two}SN_${2}.bam -normal ${2}.bam -O ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON_GR.vcf.gz --panel-of-normals ${PON} --germline-resource ${GR} --native-pair-hmm-threads 8
    # gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${3}unfiltered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON_GR.vcf.gz -O ${3}filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON_GR.vcf.gz
}

