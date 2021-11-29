#!/usr/bin/env Rscript
library(chromoMap)
library(htmlwidgets)

# chromosome files
chr_file = "length_chromosome.tsv"
# annotation files
anno_file = "D:/Hanze_Groningen/STAGE/VCF/dbSNP/R/SS6004099_merge_manual_bowtie.tsv"

head(read.table(chr_file,sep = "\t"))
head(read.table(anno_file,sep = "\t"))

chr2 = chromoMap(chr_file,anno_file,
                 data_based_color_map = T,
                 data_type = "numeric",
                 plots = "scatter", chr_width = 5, 
                 chr_length = 15,  canvas_width = 5000,
                 canvas_height = 2000)
saveWidget(widget=chr2, file="D:/Hanze_Groningen/STAGE/VCF/dbSNP/R/comp_p2.html")#, selfcontained=TRUE, knitrOptions=list())
