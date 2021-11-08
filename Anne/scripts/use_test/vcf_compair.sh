#!/usr/bin/bash

ml BCFtools/1.11-GCCcore-7.3.0

METH=bowtie
PATH_BIG=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/
# loop over samples
for sample_num in */ ; do
    if  [[ $sample_num == "S1/" ]]
    then
        echo "$sample_num"
        cd ${PATH_BIG}${sample_num}
        numbers=( )
        # loop over all folders in samples
        for d in */ ; do
            # check if folder ends with _vcf/
            if  [[ $d == *_vcf/ ]]
            then
                # name folder: 5042_vcf/
                # number = 5042
                NUM=$( echo "$(basename -- $d)" | sed 's/_vcf$//')
                # add to array
                numbers+=(${NUM})
            fi     
        done
        # Create combinations of all values in array
        declare -a encode_ids=( ${numbers[@]} )
        for (( i = 0; i < ${#encode_ids[@]}; ++i ))
        do
            for (( j = i + 1; j < ${#encode_ids[@]}; ++j ))
            do            
                PATH_DIR1=${PATH_BIG}${sample_num}/"${encode_ids[i]}"_vcf/${METH}/
                PATH_DIR2=${PATH_BIG}${sample_num}/"${encode_ids[j]}"_vcf/${METH}/
                OUTPUT_DIR=${PATH_BIG}${sample_num}/compair_"${encode_ids[i]}"_"${encode_ids[j]}"/${METH}/
                mkdir -p ${OUTPUT_DIR}
                bcftools isec ${PATH_DIR1}SS600"${encode_ids[i]}"__somatic_filtered_PON_GERM.vcf.gz ${PATH_DIR2}SS600"${encode_ids[j]}"__somatic_filtered_PON_GERM.vcf.gz -p ${OUTPUT_DIR}
            done
        done
    fi
done