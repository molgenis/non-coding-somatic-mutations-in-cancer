#!/usr/bin/bash

download_unzip_merge_files () {
    # download files: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/
    cd ./zipFile/
    wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*'
    # Loop over files
    for filename in $(ls | egrep -i 'chr[0-9,X,Y,M]{1,2}.fa.gz' )
    do 
        # Unzip files
        gunzip -c "$filename" > /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/"${filename%.*}"
        # Merge files
        cat /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/"${filename%.*}" >> /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa  
    done
    
    cd ..
    # Load packages
    ml SAMtools/1.9-foss-2018b
    ml GATK/4.1.4.1-Java-8-LTS
    # Makes dict file
    gatk CreateSequenceDictionary -R chrall.fa
    # Index file
    samtools faidx chrall.fa
}

# The file to be created in which all chromosomes are put in
FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa
# Control if file excist
if [ ! -f "$FILE" ]; then
    download_unzip_merge_files
else 
    # If the file does exist, you will be asked if you want to redo everything
    read -p "Do you want to do the download and the rest again?" yn
    case $yn in
        [Yy]* ) 
            # Remove files when 'yes' is chosen
            rm "$FILE"
            rm /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.dict
            rm /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa.fai
            rm /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/*
            rm /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/zipFile/*
            printf "DELETE SUCCESSFUL"
            printf "DOWNLOAD START"
            download_unzip_merge_files
            printf "DOWNLOAD END";;
        [Nn]* ) 
            # Does nothing when no is chosen
            printf "%s\n" "$FILE";;
        * ) 
            # When something else is chosen than yes or no, it asks to type in the right
            echo "Please answer y (=yes) or n (=no).";;
    esac
fi
