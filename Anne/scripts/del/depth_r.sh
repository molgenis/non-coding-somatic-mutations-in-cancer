#!/usr/bin/bash

# Load packages
ml SAMtools/1.9-foss-2018b
ml R/4.0.3-foss-2018b-bare


# Array of files it should download
#array=( 4109  4113  4114  4118  4119  4123  4124  4129  4128  4133  4134  4139  4138  5041  5043  5042  5044 )
array=( 4109  4113 )
# Array with file output names
#array2=( S2  S2  S2  S2  S2  S3  S3  S3  S4  S4  S4  S5  S5  S5  S6  S6  S6 )
array2=( S2  S2 )

METH=bowtie


for i in "${!array[@]}"
do
    # Number or specific tissue of a sample
    NUMBER="${array[i]}"
    # The entire file number
    FILE_NUM=SS600${NUMBER}  
    # The path where the file is located
    PATH_DIR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/"${array2[i]}"/${NUMBER}/${METH}/
    cd ${PATH_DIR}
    FILE=$(find . -type f -name "SN_S*bam")
    #https://www.biostars.org/p/104063/
    COV_FILE=${PATH_DIR}$( echo "$(basename -- $FILE)" | sed 's/.bam$//').coverage
    echo ${COV_FILE}
    samtools depth ${FILE} > ${COV_FILE}
    PATH_PLOT=${PATH_DIR}${FILE_NUM}.png
    
    # R
    echo '
    #!/usr/bin/env Rscript
    args <- commandArgs(trailingOnly=TRUE)
    print(args[1])
    print(args[2])
    coverage=read.table(args[1], sep="\t", header=F)
    print("h1")
    install.packages('reshape')
    print("h2")
    library(reshape)
    print("h3")
    coverage=rename(coverage,c(V1="Chr", V2="locus", V3="depth")) # renames the header
    print("h4")
    png(file=args[2])
    ggplot(coverage, aes(x=locus, y=depth)) +
    geom_point(colour="red", size=1, shape=20, alpha=1/3) +
    scale_y_continuous(trans = scales::log10_trans(), breaks = scales::trans_breaks("log10", function(x) 10^x))
    dev.off()
    ' > tempRscript.R
    Rscript tempRscript.R ${COV_FILE} ${PATH_PLOT}
    rm tempRscript.R

    echo "EIND"
done

