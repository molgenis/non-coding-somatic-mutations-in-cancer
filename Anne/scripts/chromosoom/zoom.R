library(karyoploteR)
library(dplyr)   

path_breast = "D:/Hanze_Groningen/STAGE/R/select_inf/vs/breast.tsv"
path_nonbreast = "D:/Hanze_Groningen/STAGE/R/select_inf/vs/nonbreast.tsv"
breast <- read.table(file = path_breast, sep = '\t', header = FALSE)
colnames(breast) <- c("gene_id", "seqnames","Start","End")
nonbreast <- read.table(file = path_nonbreast, sep = '\t', header = FALSE)
colnames(nonbreast) <- c("gene_id", "seqnames","Start","End")

# Color chromosome on genes
custom.cytobands <- toGRanges('D:/Hanze_Groningen/STAGE/db/mycytobands_R.txt')


breast$range <- paste(breast$Start, "-", breast$End)
# Make it into a GRanges, so you can make a plot with it
all_snps_breast <- GRanges(breast) #toGRanges
# Add patient ID as a column
breast$sampleName <- 'breast'

nonbreast$range <- paste(nonbreast$Start, "-", nonbreast$End)
# Make it into a GRanges, so you can make a plot with it
all_snps_nonbreast <- GRanges(nonbreast) #toGRanges
# Add patient ID as a column
nonbreast$sampleName <- 'nonbreast'

mycytobands <- read.table(file = 'D:/Hanze_Groningen/STAGE/db/mycytobands_R.txt', sep = '\t', header = TRUE)
r0 <- 0.06
r1 <- 0.50

r0_2 <- 0.56
r1_2 <- 1.00

path_info_save <- "D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/ZOOM/plotjes/"

outliers_unique_breast_path = "D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/before/outliers_unique_breast.tsv"
outliers_unique_breast <- read.table(file = outliers_unique_breast_path, sep = '\t', header = TRUE)
for(i in 1:2) { #nrow(outliers_unique_breast[:,1])
  row <- outliers_unique_breast[i,]
  print(row)
  # do stuff with row
  # Create name for the figure that saves ~12 plots      
  name <- paste(path_info_save, row$snps, '_',row$chr, '_',row$start_region, '_',row$end_region, '.png', sep="")
  # Make figure
  png(name, width = 2000, height = 1000)
  #test <- test[test$chr == 'chr13',]
  #data_interest = test[which(test$start >= 21100001)[1]:which(test$end >= 21200000)[1], ]
  data_interest = filter(mycytobands, chr == row$chr, start >= as.integer(row$start_region), end <= as.integer(row$end_region)+200000)
  
  zoom.region <- toGRanges(data.frame(row$chr, row$start_region-100000, row$end_region+100000))
  kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes=row$chr, cytobands = custom.cytobands, zoom=zoom.region)
  #y = 0.03
  #for(i in 1:nrow(data_interest)) {
  #  row_data_interest <- data_interest[i,]
  #  if (row_data_interest$start <= row$start_region-100000){
  #    kpText(kp, chr=row$chr, x=row$start_region-100000, y=y, labels=row_data_interest$name, data.panel = 1)
  #  } else if (row_data_interest$start <= row$end_region+2000) {
  #    kpText(kp, chr=row$chr, x=row$row_data_interest$start-2000, y=y, labels=row_data_interest$name, data.panel = 1)
  #  } else{
  #    kpText(kp, chr=row$chr, x=row_data_interest$start, y=y, labels=row_data_interest$name, data.panel = 1)
  #  }
  #  y = y + 0.03
  #}
  
  #r0 <- y + 0.03
  #r1 <- r0 + 0.44
  
  #  r0_2 <- r1 + 0.06
  #  r1_2 <- r0_2 + 0.44
  
  # Add base numbers: add the base numbering below each chromosome
  kpAddBaseNumbers(kp, tick.dist = 100000, minor.tick.dist = 100000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
  kp <- kpPlotDensity(kp, data=all_snps_breast, r0=r0, r1=r1,  window.size = 1000)
  # Adding Y axis to our plots may help interpreting the magnitudes of plotted data and
  # identifying the space dedicated to each plot.
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1,  cex=0.8)
  # Add labels. Labels will appear on the left side of the plot.
  kpAddLabels(kp, labels='breast', label.margin=0.04, r0=r0, r1=r1,  data.panel = 1, cex=1.2, col = "cornflowerblue")
  kp <- kpPlotDensity(kp, data=all_snps_nonbreast, r0=r0_2, r1=r1_2,  window.size = 1000)
  # Adding Y axis to our plots may help interpreting the magnitudes of plotted data and
  # identifying the space dedicated to each plot.
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=r0_2, r1=r1_2,  cex=0.8)
  # Add labels. Labels will appear on the left side of the plot.
  kpAddLabels(kp, labels='nonbreast', label.margin=0.04, r0=r0_2, r1=r1_2,  data.panel = 1, cex=1.2, col = "cornflowerblue")
  #kpAddCytobandLabels(kp, force.all = TRUE)
  kpRect(kp, chr=row$chr, x0=row$start_region, x1=row$end_region, y0=0, y1=kp$latest.plot$computed.values$max.density, col="#CCFFCBCC", data.panel="all", border=NA) #y0=0, y1=1,
  kpText(kp, chr=row$chr, x=row_data_interest$start, y=y, labels=row_data_interest$name, data.panel = 1)
  #if (nrow(data_interest) > 0){
  #  for(i in 1:nrow(data_interest)) {
  #    row_data_interest <- data_interest[i,]
  #    if (row_data_interest$start <= row$end_region+100000){
  #      print(row_data_interest$name)
  #      kpText(kp, chr=row$chr, x=row_data_interest$start, y=y, labels=row_data_interest$name, data.panel = 1)
  #    }
  #  }
  #}
  dev.off()
}












