#!/usr/bin/Rscript

library(karyoploteR)

read_file_test <- function(path){
  snps<-read.table(path,sep="\t",header=F)
  colnames(snps)<-c("gene_id", "seqnames","Start","End") 
  # Join start and end separated by a '-'
  snps$range <- paste(snps$Start, "-", snps$End)
  # Make it into a GRanges, so you can make a plot with it
  all_snps <- GRanges(snps) #toGRanges
  
  return(all_snps)
}


# Boven elkaar
kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes="chr1")
kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
data1 <- read_file_test("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/1_DO35085.bed")
kp <- kpPlotDensity(kp, data1, col="#3388FF", border="#3388FF", r0=0, r1=0.5)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5, cex=0.8)
kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='test', label.margin=0.04, r0=0, r1=0.5, data.panel = 1, cex=1, col = "cornflowerblue")

data2 <- read_file_test("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/1_DO231256.bed")
kp <- kpPlotDensity(kp, data2, col="#889F34", border="#889F34", r0=0.6, r1=1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.6, r1=1, cex=0.8)
kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0.6, r1=1)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='test2', label.margin=0.04, r0=0.6, r1=1, data.panel = 1, cex=1, col = "cornflowerblue")



# inkleuren
kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes="chr1")
kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
data1 <- read_file_test("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/1_DO35085.bed")
kp <- kpPlotDensity(kp, data1, col="#3388FF", border="#3388FF", r0=0, r1=0.5)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5, cex=0.8)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='test', label.margin=0.04, r0=0, r1=0.5, data.panel = 1, cex=1, col = "cornflowerblue")

data2 <- read_file_test("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/1_DO231256.bed")
kp <- kpPlotDensity(kp, data2, col="#889F34", border="#889F34", r0=0, r1=0.5)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5, cex=0.8)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='test2', label.margin=0.04, r0=0, r1=0.5, data.panel = 1, cex=1, col = "cornflowerblue")







# zonder inkleuren
kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes="chr1")
kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
data1 <- read_file_test("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/1_DO35085.bed")
kp <- kpPlotDensity(kp, data1, col="#FFFFFFAA", border="#3388FF", r0=0, r1=0.5)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5, cex=0.8)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='test', label.margin=0.04, r0=0, r1=0.5, data.panel = 1, cex=1, col = "cornflowerblue")

data2 <- read_file_test("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/1_DO231256.bed")
kp <- kpPlotDensity(kp, data2, col="#FFFFFFAA", border="#889F34", r0=0, r1=0.5)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5, cex=0.8)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='test2', label.margin=0.04, r0=0, r1=0.5, data.panel = 1, cex=1, col = "cornflowerblue")







########################
#
########################
kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes="chr1")
kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)

data1 <- read_file_test("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/1_DO35085.bed")
kp <- kpPlotDensity(kp, data1, col="#3388FF", border="#3388FF", r0=0, r1=0.5, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5, cex=0.8, data.panel = 1)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='test', label.margin=0.04, r0=0, r1=0.5, data.panel = 1, cex=1, col = "cornflowerblue", data.panel = 1)

data2 <- read_file_test("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/1_DO231256.bed")
kp <- kpPlotDensity(kp, data2, col="#889F34", border="#889F34", r0=0, r1=0.5, data.panel = 2)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5, cex=0.8, data.panel = 2)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='test2', label.margin=0.04, r0=0, r1=0.5, data.panel = 1, cex=1, col = "cornflowerblue", data.panel = 2)












library(TxDb.Hsapiens.UCSC.hg19.knownGene)
r0 = 0.5
r1 = 1
kp <- plotKaryotype(genome="hg19", plot.type=1, zoom="chr1:68e6-73e6")
kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
data1 <- read_file_test("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/1_DO35085.bed")
kp <- kpPlotDensity(kp, data1, col="#3388FF", border="#3388FF", r0=r0, r1=r1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1, cex=0.8)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='test', label.margin=0.04, r0=r0, r1=r1, data.panel = 1, cex=1, col = "cornflowerblue")

data2 <- read_file_test("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/1_DO231256.bed")
kp <- kpPlotDensity(kp, data2, col="#889F34", border="#889F34", r0=r0, r1=r1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1, cex=0.8)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='test2', label.margin=0.04, r0=r0, r1=r1, data.panel = 1, cex=1, col = "cornflowerblue")

genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, karyoplot = kp)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.2, cex=0.4, gene.name.position = "left")





