#!/usr/bin/bash
# Load packages
ml SAMtools/1.9-foss-2018b
cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/


# Array of files it should download
array=( ./SS6004099.bam  ./SS6004094.bam  ./SS6004104.bam  ./SS6004113.bam  ./SS6004114.bam
./SS6004118.bam  ./SS6004109.bam  ./SS6004119.bam  ./SS6004124.bam  ./SS6004123.bam
./SS6004129.bam  ./SS6004133.bam  ./SS6004128.bam  ./SS6004134.bam  ./SS6004139.bam
./SS6004138.bam  ./SS6005041.bam  ./SS6005043.bam  ./SS6005042.bam  ./SS6005044.bam )
# Array with file output names
array2=( ./SS6004099.sorted.bam  ./SS6004094.sorted.bam  ./SS6004104.sorted.bam  ./SS6004113.sorted.bam  
        ./SS6004114.sorted.bam ./SS6004118.sorted.bam  ./SS6004109.sorted.bam  ./SS6004119.sorted.bam  
        ./SS6004124.sorted.bam  ./SS6004123.sorted.bam ./SS6004129.sorted.bam  ./SS6004133.sorted.bam
        ./SS6004128.sorted.bam  ./SS6004134.sorted.bam  ./SS6004139.sorted.bam ./SS6004138.sorted.bam  
        ./SS6005041.sorted.bam  ./SS6005043.sorted.bam  ./SS6005042.sorted.bam  ./SS6005044.sorted.bam )

for i in "${!array[@]}"
do
    # Checked if the sorted BAM file already exists.
    while [ ! -f "${array2[i]}" ]
    do
        # If this file does not yet exist, the original BAM file will be sorted and indexed.
        printf "sort: %s" "${array[i]}"
        samtools sort "${array[i]}" -o "${array2[i]}"
        printf "index: %s" "${array2[i]}"
        samtools index "${array2[i]}"
        printf "Done with: %s" "${array[i]}"
    done
done