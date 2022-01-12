#!/usr/bin/bash

# Load package
ml BCFtools/1.11-GCCcore-7.3.0
ml GATK/4.1.4.1-Java-8-LTS

# Which method
METH=bowtie
# The general path
PATH_BIG=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/
# Loop over samples (folders S1, S2 etc.)
for sample_num in */ ; do
    if  [[ $sample_num == "S1/" ]]
    then
        echo "$sample_num"
        cd ${PATH_BIG}${sample_num}
        # Empty array
        numbers=( )
        # Loop over all folders in sample
        for d in */ ; do
            # Check if folder ends with _vcf/
            if  [[ $d == *_vcf/ ]]
            then
                # Name folder: 5042_vcf/
                # NUM = 5042
                NUM=$( echo "$(basename -- $d)" | sed 's/_vcf$//')
                # Add to array
                numbers+=(${NUM})
            fi     
        done
        # Create combinations of all values in array
        declare -a ids=( ${numbers[@]} )
        for (( i = 0; i < ${#ids[@]}; ++i ))
        do
            for (( j = i + 1; j < ${#ids[@]}; ++j ))
            do            
                # Path for first number
                PATH_DIR1=${PATH_BIG}${sample_num}/"${ids[i]}"_vcf/${METH}/
                # Path for second number
                PATH_DIR2=${PATH_BIG}${sample_num}/"${ids[j]}"_vcf/${METH}/
                # Output path
                OUTPUT_DIR=${PATH_BIG}${sample_num}/compair_"${ids[i]}"_"${ids[j]}"/${METH}/
                # Make folder
                mkdir -p ${OUTPUT_DIR}
                # Compare two VCF files
                bcftools isec ${PATH_DIR1}SS600"${ids[i]}"__somatic_filtered_PON_GERM.vcf.gz ${PATH_DIR2}SS600"${ids[j]}"__somatic_filtered_PON_GERM.vcf.gz -p ${OUTPUT_DIR}
                for filename in ${OUTPUT_DIR}*.vcf; do
                    bcftools view ${filename} -Oz -o ${filename}.gz
                    bcftools index ${filename}.gz
                #     gatk IndexFeatureFile -I ${filename}
                done
            done
        done
    fi
done