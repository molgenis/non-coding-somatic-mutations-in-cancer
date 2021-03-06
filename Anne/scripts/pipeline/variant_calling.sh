#!/usr/bin/bash

# bam --> VCF
mutect2_vcf() {
    # Mutect2:            Mutect2 is designed to call somatic variants only
    # FilterMutectCalls:  Filter somatic SNVs and indels called by Mutect2
    # Do this now in 4 different ways, so that you can check later whether there is much difference.
    mkdir -p ${1}
    # sample name
    echo "####sample name\n"
    gatk Mutect2 -R ${PATH_GENOOM}${CHROM}.fa ${2} -O ${3}somatic_unfiltered.vcf.gz --native-pair-hmm-threads 8
    gatk FilterMutectCalls -R ${PATH_GENOOM}${CHROM}.fa -V ${3}somatic_unfiltered.vcf.gz -O ${3}somatic_filtered.vcf.gz
    # PON
    #echo "####PON\n"
    #gatk Mutect2 -R ${PATH_GENOOM}${CHROM}.fa ${2} -O ${3}somatic_unfiltered_PON.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8
    #gatk FilterMutectCalls -R ${PATH_GENOOM}${CHROM}.fa -V ${3}somatic_unfiltered_PON.vcf.gz -O ${3}somatic_filtered_PON.vcf.gz
    #GR
    #echo "####GR\n"
    #gatk Mutect2 -R ${PATH_GENOOM}${CHROM}.fa ${2} -O ${3}somatic_unfiltered_GERM.vcf.gz --native-pair-hmm-threads 8 --germline-resource ${GR}
    #gatk FilterMutectCalls -R ${PATH_GENOOM}${CHROM}.fa --V ${3}somatic_unfiltered_GERM.vcf.gz -O ${3}somatic_filtered_GERM.vcf.gz
    #PON and GR
    echo "####PON and GR\n"
    gatk Mutect2 -R ${PATH_GENOOM}${CHROM}.fa ${2} -O ${3}somatic_unfiltered_PON_GERM.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8 --germline-resource ${GR}
    gatk FilterMutectCalls -R ${PATH_GENOOM}${CHROM}.fa --V ${3}somatic_unfiltered_PON_GERM.vcf.gz -O ${3}somatic_filtered_PON_GERM.vcf.gz
}

