#!/usr/bin/bash

#sets working directory
#cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/

# array of files it should download
array=( EGAF00000151605  EGAF00000151604  EGAF00000151606  EGAF00000151608  EGAF00000151609
EGAF00000151610  EGAF00000151607  EGAF00000151611  EGAF00000151613  EGAF00000151612
EGAF00000151614  EGAF00000151616  EGAF00000151615  EGAF00000151617  EGAF00000151619
EGAF00000151618  EGAF00000151620  EGAF00000151622  EGAF00000151621  EGAF00000151623 )
# array with file output names
array2=( ./SS6004099.bam  ./SS6004094.bam  ./SS6004104.bam  ./SS6004113.bam  ./SS6004114.bam
./SS6004118.bam  ./SS6004109.bam  ./SS6004119.bam  ./SS6004124.bam  ./SS6004123.bam
./SS6004129.bam  ./SS6004133.bam  ./SS6004128.bam  ./SS6004134.bam  ./SS6004139.bam
./SS6004138.bam  ./SS6005041.bam  ./SS6005043.bam  ./SS6005042.bam  ./SS6005044.bam )

for i in "${!array[@]}"
do
    #printf "%s is in %s\n" "${array[i]}" "${array2[i]}"
    while [ ! -f "${array2[i]}" ]
    do
        printf "HALLOOOOOOOOOOOO %s is in %s\n" "${array[i]}" "${array2[i]}i"
        #https://github.com/EGA-archive/ega-download-client
        pyega3 -c 30 -cf ~/ega_credentials.json fetch "${array[i]}" --saveto "${array2[i]}"
        #touch "${array2[i]}"
    done
done




