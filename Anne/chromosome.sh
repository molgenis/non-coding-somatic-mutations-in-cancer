#!/usr/bin/bash

download_unzip_merge_files () {
    # download files: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/
    cd ./zipFile/
    wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*'
    # Loop over files
    for filename in $(ls | egrep -i 'chr[0-9,X,Y,M]{1,2}.fa.gz' )
    do 
        #printf "%s\n" "$filename"
        # unzip files
        gunzip -c "$filename" > /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/"${filename%.*}"
        # merge files
        cat /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/"${filename%.*}" >> /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa  
    done
    
    cd ..
    ml SAMtools/1.9-foss-2018b
    ml GATK/4.1.4.1-Java-8-LTS
    # dict file
    gatk CreateSequenceDictionary -R chrall.fa
    # index file
    samtools faidx chrall.fa
}


FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa
# Control if file excist
if [ ! -f "$FILE" ]; then
    download_unzip_merge_files
else 
    read -p "Do you want to do the download and the rest again?" yn
    case $yn in
        [Yy]* ) 
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
            printf "%s\n" "$FILE";;
        * ) 
            echo "Please answer y (=yes) or n (=no).";;
    esac
fi