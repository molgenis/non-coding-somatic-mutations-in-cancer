#https://bernatgel.github.io/karyoploter_tutorial//Examples/GeneDensityIdeograms/GeneDensityIdeograms.html
#https://bernatgel.github.io/karyoploter_tutorial/Tutorial/PlotDensity/PlotDensity.html
library(karyoploteR)
# custom genomes
#custom.genome <- toGRanges(data.frame(chr=c("A", "B"), start=c(1, 1), end=c(100, 200)))
#kp <- plotKaryotype(genome = custom.genome)

#custom.genome <- toGRanges("Tutorial/CustomGenomes/mygenome.txt")
#kp <- plotKaryotype(genome = custom.genome)

#custom.genome <- toGRanges("Tutorial/CustomGenomes/mygenome.txt")
#custom.cytobands <- toGRanges("Tutorial/CustomGenomes/mycytobands.txt")
#kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands)


create_data_karyoploteR <- function(file_name) {
  path = 'D:/Hanze_Groningen/STAGE/bed/files/'
  snps22<-read.table(paste(path, file_name, sep=""),sep="\t",header=F)
  colnames(snps22)<-c("gene_id", "seqnames","Start","End")
  snps22$range <- paste(snps22$Start, "-", snps22$End)
  genes <- GRanges(snps22) #toGRanges
  return(genes)
}

add_density <- function(kp, data, r0, r1, label) {
  kp <- kpPlotDensity(kp, data=data, r0=r0, r1=r1, window.size = 0.5e6) #window.size = 0.5e6 100000
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1, cex=0.8)
  kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1)
  kpAddLabels(kp, labels=label, label.margin=0.04, r0=r0, r1=r1, data.panel = 1, cex=0.8, col = "cornflowerblue")
  return(kp)
}


kp <- plotKaryotype()
# select chromomes
kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chr10", "chr12", "chr22"), main="Gene Density")
# add base numbers
kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")
# add Density

sample1 <- create_data_karyoploteR('SS6004099_merge_manual_bowtie.tsv')
kp <- add_density(kp, sample1, 0, 0.27, 'SS6004099')

sample2 <- create_data_karyoploteR('SS6004104_merge_manual_bowtie.tsv')
kp <- add_density(kp, sample2, 0.34, 0.61, 'SS6004104')

sample3 <- create_data_karyoploteR('SS6004113_merge_manual_bowtie.tsv')
kp <- add_density(kp, sample3, 0.68, 0.95, 'SS6004113')

sample4 <- create_data_karyoploteR('SS6004114_merge_manual_bowtie.tsv')
kp <- add_density(kp, sample4, 1.02, 1.29, 'SS6004114')
