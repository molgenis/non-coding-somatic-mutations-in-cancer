#!/usr/bin/bash

FILE=bowtie_S1_numT_1_numHC_1_SS6004094_SS6004099.txt
INPUT_FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/S1/${FILE}
arg1=$(sed '1q;d' ${INPUT_FILE})
arg2=$(sed '2q;d' ${INPUT_FILE})
arg3=$(sed '3q;d' ${INPUT_FILE})

source /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/make_vcf_auto.sh

mutect2_vcf ${arg1} "${arg2}" ${arg3}
