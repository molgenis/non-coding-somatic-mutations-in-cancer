#!/usr/bin/bash

#sets working directory
#cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/

# array of files it should download
array=( EGAF00000151604 EGAF00000151605 EGAF00000151606 )
# array with file output names
array2=( ./SS6004094.bam ./SS6004099.bam ./SS6004104.bam )

for i in "${!array[@]}"
do
    #printf "%s is in %s\n" "${array[i]}" "${array2[i]}"
    while [ ! -f "${array2[i]}" ]
    do
        #printf "%s is in %s\n" "${array[i]}" "${array2[i]}i"
        pyega3 -cf ~/ega_credentials.json fetch "${array[i]}" --saveto "${array2[i]}" 
        #touch "${array2[i]}"
    done
done




