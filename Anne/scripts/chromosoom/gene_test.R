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
all_snps <- GRanges(breast) #toGRanges
# Add patient ID as a column
breast$sampleName <- 'breast'


test <- read.table(file = 'D:/Hanze_Groningen/STAGE/db/mycytobands_R.txt', sep = '\t', header = TRUE)
test <- test[test$chr == 'chr17',]
data_interest = test[which(test$start >= 41196311)[1]:which(test$end >= 41277500)[1], ]
print(data_interest)


zoom.region <- toGRanges(data.frame('chr17', 41096311, 41377500))
kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes='chr17', cytobands = custom.cytobands, zoom=zoom.region)
kp <- kpPlotDensity(kp, data=all_snps, window.size = 1000)
# Add base numbers: add the base numbering below each chromosome
kpAddBaseNumbers(kp, tick.dist = 100000, minor.tick.dist = 100000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
# Adding Y axis to our plots may help interpreting the magnitudes of plotted data and
# identifying the space dedicated to each plot.
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.8)
# Add labels. Labels will appear on the left side of the plot.
kpAddLabels(kp, labels='label', label.margin=0.04, data.panel = 1, cex=1.2, col = "cornflowerblue")
#kpAddCytobandLabels(kp, force.all = TRUE)
#kpRect(kp, chr="chr13", x0=21100001, x1=21200000, y0=0, y1=kp$latest.plot$computed.values$max.density, col="#AAFFCBDD", data.panel="all", border=NA) #y0=0, y1=1, 
kpText(kp, chr="chr17", x=data_interest$start, y=0.5, labels=data_interest$name, data.panel = 1)

