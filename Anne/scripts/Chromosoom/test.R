library(karyoploteR)

create_data_karyoploteR <- function(path) {
  snps22<-read.table(path,sep="\t",header=F)
  colnames(snps22)<-c("gene_id", "seqnames","Start","End")
  snps22$range <- paste(snps22$Start, "-", snps22$End)
  genes <- GRanges(snps22) #toGRanges
  return(genes)
}

add_density <- function(kp, data, r0, r1, label) {
  kp <- kpPlotDensity(kp, data=data, r0=r0, r1=r1, window.size = 100000) #window.size = 0.5e6 100000
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1, cex=0.8)
  #kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1)
  kpAddLabels(kp, labels=label, label.margin=0.04, r0=r0, r1=r1, data.panel = 1, cex=1.2, col = "cornflowerblue")
  return(kp)
}

numbers <-  seq(1, 22, by=1)
chrom <- append(numbers, c('X', 'Y'))

filenames <- list.files("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/", pattern="*.bed", full.names=TRUE)
for (chr in chrom){
  chr <- paste('chr', chr, sep="")
  for (i in 1:length(filenames)){
      r0=0
      r1=0.27

    if ( (i%%5) == 0 || i == 1){
      if (i != 1){
        dev.off()
      }
      counter <- 0
      basename <- strsplit(basename(filenames[i]), "_", fixed = TRUE)[[1]][1]
      name <- paste('D:/Hanze_Groningen/STAGE/bed/PLOTS/kary/', chr, '/plot', i, basename, '.png', sep="")
      png(name, width = 2000, height = 1000) #width = 2000, height = 1000
      kp <- plotKaryotype()
      # select chromomes
      kp <- plotKaryotype(genome="hg19", plot.type=1, main="Gene Density", chromosomes=c(chr)) #chromosomes=c("chr22")
      # add base numbers
      kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
      sample1 <- create_data_karyoploteR(filenames[i])
      add <- (counter*0.30)
      kp <- add_density(kp, sample1, r0+add, r1+add, basename)

    } else{
      basename <- strsplit(basename(filenames[i]), "_", fixed = TRUE)[[1]][1]
      sample1 <- create_data_karyoploteR(filenames[i])
      counter <- counter + 1
      add <- (counter*0.30)
      kp <- add_density(kp, sample1, r0+add, r1+add, basename)

    }
    if (i == length(filenames)) {
      dev.off()
    }
  }
}


