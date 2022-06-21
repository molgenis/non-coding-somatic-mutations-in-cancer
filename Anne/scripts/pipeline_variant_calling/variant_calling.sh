#!/usr/bin/bash

# Function that runs Mutect2
# Mutect2 is run first without additional parameters.
# Then Mutect2 is run with PoN and Germline resources
# A vcf file comes out with somatic mutations of the sample

mutect2_vcf() {
    ml GATK/4.1.4.1-Java-8-LTS
    # Mutect2:            Mutect2 is designed to call somatic variants only
    # FilterMutectCalls:  Filter somatic SNVs and indels called by Mutect2
    # Do this now in 2 different ways, so that you can check later whether there is much difference.
    mkdir -p ${1}
    # sample name
    echo "####sample name\n"
    gatk Mutect2 -R ${PATH_GENOOM}${CHROM}.fa ${2} -O ${3}somatic_unfiltered.vcf.gz --native-pair-hmm-threads 8
    gatk FilterMutectCalls -R ${PATH_GENOOM}${CHROM}.fa -V ${3}somatic_unfiltered.vcf.gz -O ${3}somatic_filtered.vcf.gz
    #PON and GR
    echo "####PON and GR\n"
    gatk Mutect2 -R ${PATH_GENOOM}${CHROM}.fa ${2} -O ${3}somatic_unfiltered_PON_GERM.vcf.gz --panel-of-normals ${PON} --native-pair-hmm-threads 8 --germline-resource ${GR}
    gatk FilterMutectCalls -R ${PATH_GENOOM}${CHROM}.fa --V ${3}somatic_unfiltered_PON_GERM.vcf.gz -O ${3}somatic_filtered_PON_GERM.vcf.gz
}

