#"D:/Hanze_Groningen/STAGE/bed/snp132_ucsc_hg19.bed"

#name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID	alignID
library(data.table) #setnames
library("chromPlot") #Plot chromosoom

delete_chr <- function(data) {
  # Delete 'chr' for each chromosome
  data$chr <- gsub('chr', '', data$chr)
  #merge_gene_exon$chr <- sapply(strsplit(merge_gene_exon$chr,"_"), `[`, 1)
  # Replace some letters (X, Y, M) with numbers
  data$chr<-replace(data$chr, data$chr=='X', 23)
  data$chr<-replace(data$chr, data$chr=='Y', 24)
  data$chr<-replace(data$chr, data$chr=='M', 25)
  # Make integers
  data$chr <- as.integer(data$chr)
  data$Start <- as.integer(data$Start)
  data$End <- as.integer(data$End)
  # Order on columns
  sort_merge <- data[with(data, order(chr, Start, End)),]
  #Remove NA values
  delete_na_values <- sort_merge[!is.na(sort_merge$chr), ]
  delete_na_values$Chrom <- delete_na_values$chr
  # Replace the numbers with the letters again
  delete_na_values$Chrom<-replace(delete_na_values$Chrom, delete_na_values$Chrom==23, 'X')
  delete_na_values$Chrom<-replace(delete_na_values$Chrom, delete_na_values$Chrom==24, 'Y')
  delete_na_values$Chrom<-replace(delete_na_values$Chrom, delete_na_values$Chrom==25, 'M')
  
  return(delete_na_values)
}

##########
# gaps
# https://bioconductor.riken.jp/packages/3.6/bioc/vignettes/chromPlot/inst/doc/chromPlot.pdf
# The format for the ‘Gap’track in the Table Browser of the UCSC website
# From this table, chromPlot extracts the number of chromosomes, chromosomes names and lengths, and the
# position of centromeres (shown as solid circles).
##########
# Zelfde als hg_gap?
#chromosomes<-read.table("D:/Hanze_Groningen/STAGE/bed/snp132_ucsc_hg19_gap.bed",sep="\t",header=F)
# bin	chrom	chromStart	chromEnd	ix	n	size	type	bridge
#colnames(chromosomes)<-c("bin", "chr","Start", "End", 'ix',	'n',	'size',	'Name',	'bridge')
#chromosomes <- chromosomes[,c("chr","Start", "End", "Name")]
#sort_chromosomes <- delete_chr(chromosomes)
data(hg_gap)
head(hg_gap)


##########
# bands
# For the colors in the chromosome
##########
# Open and read files
# exon file was edited with python script

gene_exon_files <- function(file_name) {
  # Remove '#' in file, so that the file has headers
  path <- "D:/Hanze_Groningen/STAGE/bed/"
  exon <-read.table(paste(path, file_name, '_exon.bed', sep=""),sep="\t",header=T)
  gene <-read.table(paste(path, file_name, '_gene.bed', sep=""),sep="\t",header=T)
  # Select columns
  filter_exon <- exon[,c("name","chrom","strand", "exonStarts", "exonEnds")]
  filter_gene <- gene[,c("name","chrom","strand", "txStart", "txEnd")]
  # Rename columns
  setnames(filter_exon, old = c('exonStarts','exonEnds', 'chrom'), new = c('Start','End', 'chr'))
  setnames(filter_gene, old = c('txStart','txEnd', 'chrom'), new = c('Start','End', 'chr'))
  # Add column Group
  filter_exon$Group <- 'exon'
  filter_gene$Group <- 'gene'
  # Merge files
  merge_gene_exon <- rbind(filter_gene,filter_exon)
  sort_gene_exon <- delete_chr(merge_gene_exon)
  # Make separate values of genes and exons again
  only_gene <- sort_gene_exon[sort_gene_exon$Group == 'gene',]
  only_exon_in_gene <- sort_gene_exon[sort_gene_exon$Group == 'exon',]
  return(only_gene)
}

UCSC_gene <- gene_exon_files('snp132_ucsc_hg19_UCSC')
NCBI_gene <- gene_exon_files('snp132_ucsc_hg19_NCBI')
V38_gene <- gene_exon_files('snp132_ucsc_hg19_V38')



##########
# annot1
# For the histogram
##########
# The 'somatic' mutations out the vcf files

create_data_chromplot <- function(file_name) {
  path = 'D:/Hanze_Groningen/STAGE/bed/files/'
  snps<-read.table(paste(path, file_name, sep=""),sep="\t",header=F)
  colnames(snps)<-c("Name", "chr","Start","End")
  sort_snps <- delete_chr(snps)
  return(sort_snps)
}
chro_sample1 <- create_data_chromplot('SS6004099_merge_manual_bowtie.tsv')

chro_sample2 <- create_data_chromplot('SS6004104_merge_manual_bowtie.tsv')

chro_sample3 <- create_data_chromplot('SS6004113_merge_manual_bowtie.tsv')

chro_sample4 <- create_data_chromplot('SS6004114_merge_manual_bowtie.tsv')




##########
# ChromPlot
# Making the chromosome plots
##########
#chromPlot(gaps=delete_na_values, bands=only_gene, figCols=6)  
chromPlot(gaps=hg_gap, bands=UCSC_gene, colBand = "##9e9f93", annot1=chro_sample1,  colAnnot1="brown", figCols=3, scale.title="Counts", title = "UCSC_gene", legChrom=22, chr=c(20:22))

#chromPlot(gaps=hg_gap, bands=NCBI_gene, colBand = "##9e9f93", annot1=chro_sample1,  colAnnot1="brown", figCols=3, scale.title="Counts", title = "NCBI_gene", legChrom=22, chr=c(20:22))

#chromPlot(gaps=hg_gap, bands=V38_gene, colBand = "##9e9f93", annot1=chro_sample1,  colAnnot1="brown", figCols=3, scale.title="Counts", title = "V38_gene", legChrom=22, chr=c(20:22))

chromPlot(gaps=hg_gap, bands=UCSC_gene, colBand = "##9e9f93", annot1=chro_sample2,  colAnnot1="brown", figCols=3, scale.title="Counts", title = "UCSC_gene", legChrom=22, chr=c(20:22))

chromPlot(gaps=hg_gap, bands=UCSC_gene, colBand = "##9e9f93", annot1=chro_sample3,  colAnnot1="brown", figCols=3, scale.title="Counts", title = "UCSC_gene", legChrom=22, chr=c(20:22))

chromPlot(gaps=hg_gap, bands=UCSC_gene, colBand = "##9e9f93", annot1=chro_sample4,  colAnnot1="brown", figCols=3, scale.title="Counts", title = "UCSC_gene", legChrom=22, chr=c(20:22))