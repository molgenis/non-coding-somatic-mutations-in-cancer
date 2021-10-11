#!/usr/bin/bash

cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals
ml BCFtools/1.11-GCCcore-7.3.0

#FILE="somatic-b37_Mutect2-WGS-panel-b37.vcf"
FILE="somatic-hg38_1000g_pon.hg38.vcf"

# https://www.biostars.org/p/160863/
bcftools view --header-only ${FILE} | sed 's/##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/' | sed 's/##contig=<ID=MT/##contig=<ID=chrM/' > header_${FILE}
bcftools view --no-header ${FILE} | sed 's/^/chr/' | sed 's/^MT/chrM/' >> noHeader_${FILE}
cat header_${FILE} noHeader_${FILE} >> merge_${FILE}