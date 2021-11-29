library(ggplot2) 

path1 = "D:/Hanze_Groningen/STAGE/VCF/dbSNP/R_SS6004099_merge_manual_bowtie.tsv"

dataframe <- read.table(path1,sep = "\t", header=TRUE)
dataframe <- dataframe[,-(4:5)]
head(dataframe)
ggplot(dataframe) + geom_histogram(aes(x=Start),binwidth=1e6) + facet_grid(Chr ~Sample)