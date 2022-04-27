library(karyoploteR)

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


test <- read.table(file = 'D:/Hanze_Groningen/STAGE/db/mycytobands_R.txt', sep = '\t', header = TRUE)
test <- test[test$chr == 'chr13',]
data_interest = test[which(test$start >= 21100001)[1]:which(test$end >= 21200000)[1], ]
print(data_interest)

r0 <- 0.06
r1 <- 0.50

r0_2 <- 0.56
r1_2 <- 1.00

zoom.region <- toGRanges(data.frame('chr13', 21000001, 21300000))
kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes='chr13', cytobands = custom.cytobands, zoom=zoom.region)
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
kpRect(kp, chr="chr13", x0=21100001, x1=21200000, y0=0, y1=kp$latest.plot$computed.values$max.density, col="#AAFFCBDD", data.panel="all", border=NA) #y0=0, y1=1, 
kpText(kp, chr="chr13", x=data_interest$start+10000, y=0.03, labels=data_interest$name, data.panel = 1)

