#"D:/Hanze_Groningen/STAGE/bed/snp132_ucsc_hg19.bed"
# remove # from file headers #name
#name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID	alignID
library(data.table) #setnames
library("chromPlot") #Plot chromosoom

exon <-read.table("D:/Hanze_Groningen/STAGE/bed/snp132_ucsc_hg19_exon.bed",sep="\t",header=T)
gene <-read.table("D:/Hanze_Groningen/STAGE/bed/snp132_ucsc_hg19_gene.bed",sep="\t",header=T)

filter_exon <- exon[,c("name","chrom","strand", "exonStarts", "exonEnds", "proteinID", "alignID")]
setnames(filter_exon, old = c('exonStarts','exonEnds'), new = c('Start','End'))
filter_exon$Group <- 'exon'

filter_gene <- gene[,c("name","chrom","strand", "txStart", "txEnd", "proteinID", "alignID")]
setnames(filter_gene, old = c('txStart','txEnd'), new = c('Start','End'))
filter_gene$Group <- 'gene'
#cancate
merge_gene_exon <- rbind(filter_gene,filter_exon)
merge_gene_exon$chr <- gsub('chr', '', merge_gene_exon$chrom)
#merge_gene_exon$chr <- sapply(strsplit(merge_gene_exon$chr,"_"), `[`, 1)
#replaces the negative numbers with zeros
merge_gene_exon$chr<-replace(merge_gene_exon$chr, merge_gene_exon$chr=='X', 23)
merge_gene_exon$chr<-replace(merge_gene_exon$chr, merge_gene_exon$chr=='Y', 24)
merge_gene_exon$chr<-replace(merge_gene_exon$chr, merge_gene_exon$chr=='M', 25)
#merge_gene_exon$chr<-replace(merge_gene_exon$chr, merge_gene_exon$chr=='Un', 26)
#Make integers 
merge_gene_exon$chr <- as.integer(merge_gene_exon$chr)
merge_gene_exon$Start <- as.integer(merge_gene_exon$Start)
merge_gene_exon$End <- as.integer(merge_gene_exon$End)
#Order on columns
sort_merge <- merge_gene_exon[with(merge_gene_exon, order(chr, Start, End)),]
#Remove NA values
delete_na_values <- sort_merge[!is.na(sort_merge$chr), ]
#setnames(delete_na_values, old = c('chrom', 'name'), new = c('Chrom', 'Name'))
delete_na_values$Chrom <- delete_na_values$chr
delete_na_values$Chrom<-replace(delete_na_values$Chrom, delete_na_values$Chrom==23, 'X')
delete_na_values$Chrom<-replace(delete_na_values$Chrom, delete_na_values$Chrom==24, 'Y')
delete_na_values$Chrom<-replace(delete_na_values$Chrom, delete_na_values$Chrom==25, 'M')

only_gene <- delete_na_values[delete_na_values$Group == 'gene',]
only_gene$Color <- 'black'
only_exon_in_gene <- delete_na_values[delete_na_values$Group == 'exon',]

snps<-read.table("SS6004094_merge_mutect2_bowtie.tsv",sep="\t",header=F)
colnames(snps)<-c("Name", "Chrom","Start","End")


#chromPlot(gaps=delete_na_values, bands=only_gene, figCols=6)  
chromPlot(gaps=delete_na_values, bands=only_gene, colBand = "##9e9f93", annot1=snps,  colAnnot1="brown", chr=c(20:22), figCols=3, scale.title="Counts", title = "Chromosome Plot", legChrom=22)















