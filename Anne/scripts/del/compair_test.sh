#!/usr/bin/bash

ml BCFtools/1.11-GCCcore-7.3.0

PATH_BIG=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/
# loop over samples
for sample_num in */ ; do
    echo "$sample_num"
    cd ${PATH_BIG}${sample_num}
    numbers=( )
    # loop over all folders in samples
    for d in */ ; do
        # check if folder ends with _vcf/
        if  [[ $d == *_vcf/ ]]
        then
            # add to array
            numbers+=("${d}")
        fi     
    done
    # Loop over array
    echo ${numbers[@]}
    declare -a encode_ids=( ${numbers[@]} )
    for (( i = 0; i < ${#encode_ids[@]}; ++i )); do
        for (( j = i + 1; j < ${#encode_ids[@]}; ++j )); do
          echo "${encode_ids[i]}"
          echo "${encode_ids[j]}"
          echo "----"
        done
    done
done