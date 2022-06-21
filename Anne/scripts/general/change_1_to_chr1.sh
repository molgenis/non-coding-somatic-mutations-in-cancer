#!/usr/bin/bash

cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals
# Load packages
ml BCFtools/1.11-GCCcore-7.3.0

FILE="somatic-b37_Mutect2-WGS-panel-b37.vcf"

# https://www.biostars.org/p/160863/
# Convert all chromosomes that are only called by a number to chr with the number so for example 1 now becomes chr1
# First this is done in the header.
# Then this is done in everything except the header.
bcftools view --header-only ${FILE} | sed 's/##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/' | sed 's/##contig=<ID=MT/##contig=<ID=chrM/' > header_${FILE}
bcftools view --no-header ${FILE} | sed 's/^/chr/' | sed 's/^MT/chrM/' >> noHeader_${FILE}
cat header_${FILE} noHeader_${FILE} >> merge_${FILE}