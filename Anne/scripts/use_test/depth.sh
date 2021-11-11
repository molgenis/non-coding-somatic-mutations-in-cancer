#!/usr/bin/bash

# Load packages
ml SAMtools/1.9-foss-2018b
ml R/4.0.3-foss-2018b-bare


# Array of files it should download
#array=( 4094  4099  4104  4109  4113  4114  4118  4119  4123  4124  4129  4128  4133  4134  4139  4138  5041  5043  5042  5044 )
array=( 4109   )
# Array with file output names
#array2=( S1  S1  S1  S2  S2  S2  S2  S2  S3  S3  S3  S4  S4  S4  S5  S5  S5  S6  S6  S6 )
array2=( S2  )
PATH_GLOBAL=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/
OUTPUT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/depth/
METH=bowtie


for i in "${!array[@]}"
do
    # Number or specific tissue of a sample
    NUMBER="${array[i]}"
    # The entire file number
    FILE_NUM=SS600${NUMBER}  
    # The path where the file is located
    PATH_DIR=${PATH_GLOBAL}"${array2[i]}"/${NUMBER}/${METH}/
    cd ${PATH_DIR}
    FILE=$(find . -type f -name "SN_S*bam")
    #https://www.biostars.org/p/104063/
    COV_FILE=${PATH_DIR}$( echo "$(basename -- $FILE)" | sed 's/.bam$//').coverage
    echo ${COV_FILE}
    samtools depth ${FILE} > ${COV_FILE}
    PATH_PLOT=${OUTPUT}"${array2[i]}"_${FILE_NUM}_${METH}.png
    
    
    python3 ${PATH_GLOBAL}depth_plot.py ${COV_FILE} ${PATH_PLOT}

    echo "EIND"
done

