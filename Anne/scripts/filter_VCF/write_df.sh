#!/usr/bin/bash

# Load packages
ml BCFtools/1.11-GCCcore-7.3.0
ml Anaconda3/5.3.0
source activate stage

PATH_PLOT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/FORMAT/

> /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/FORMAT/part_command.txt

# Loop over lines in file
while IFS= read -r line; do
  # Call python file
  # The python file creates dataframe from the manually created tumor files 
  # (by comparing hc against tumor, and grabbing the unique tumor SNPs)
  python3 ${PATH_PLOT}make_df.py ${line} ${PATH_PLOT} ${PATH_PLOT}part_command.txt
done < "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/files.txt"

echo 'LOOP OVER COMMANDS'
while IFS= read -r line; do
  OUTPUT_FILE=$(echo $line | rev | cut -d ' ' -f '1' | rev)
  DIR=${OUTPUT_FILE%/*}
  mkdir -p ${DIR}
  #echo "$(basename -- $FILE_NAME)"
  echo $line
  bcftools query -Hf $line
done < "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/FORMAT/part_command.txt"

echo 'END'