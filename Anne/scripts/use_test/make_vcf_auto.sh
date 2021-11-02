#!/usr/bin/bash

# Path to the reference genome.
PATH_REF=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/
# File name of the reference genome
REF=chr22.fa
# Path to the file of the panel of normals (PoN)
PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/merge_somatic-b37_Mutect2-WGS-panel-b37.vcf
# Path to the file of the germline resource (GR)
GR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GermlineResource/merge_af-only-gnomad.raw.sites.vcf

# load packages
ml Anaconda3/5.3.0
ml SAMtools/1.9-foss-2018b
ml GATK/4.1.4.1-Java-8-LTS
ml BCFtools/1.11-GCCcore-7.3.0
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS

# bam --> VCF
mutect2_vcf() {
    '''
    Mutect2:            Mutect2 is designed to call somatic variants only
    FilterMutectCalls:  Filter somatic SNVs and indels called by Mutect2
    Do this now in 4 different ways, so that you can check later whether there is much difference.
    '''
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

