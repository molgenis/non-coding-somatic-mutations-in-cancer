#!/usr/bin/env Rscript
library("chromPlot")

snps<-read.table("SS6004094_merge_mutect2_bowtie.tsv",sep="\t",header=F)
colnames(snps)<-c("Name", "Chrom","Start","End")

snps2<-read.table("SS6004094_merge_mutect2_bowtie.tsv",sep="\t",header=F, stringsAsFactors=FALSE)
colnames(snps2)<-c("Name", "Chrom","Start","End")

# put the chromosomes in the good order: chr1, chr2, chr22, chrX
#goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")
#snps$chr <- factor(snps$chr,levels=goodChrOrder)

#https://bioconductor.riken.jp/packages/3.6/bioc/vignettes/chromPlot/inst/doc/chromPlot.pdf

#chromPlot(gaps=snps, annot1=snps, chr=c(22))
#TODO uitzoeken waar deze vandaan komen
data(hg_cytoBandIdeo)
head(hg_cytoBandIdeo)
#chromPlot(bands=hg_cytoBandIdeo, gaps=hg_gap, chr=c(22), figCols=6)
#snps2$Colors <- 'red'
#chromPlot(gaps=snps, bands=snps2, chr=c(22), figCols=3)

chromPlot(gaps=snps, bands=hg_cytoBandIdeo, annot1=snps, chr=c(22), figCols=3, scale.title="Counts", legChrom='22', title = "Chromosome Plot")

chromPlot(gaps=snps, annot1=snps, legChrom=22)



