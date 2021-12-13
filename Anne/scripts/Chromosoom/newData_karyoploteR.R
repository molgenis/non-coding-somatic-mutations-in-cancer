library(karyoploteR)
#library(zeallot) #https://stackoverflow.com/questions/13538628/returning-different-data-frames-in-a-function-r

create_data_karyoploteR <- function(path, df, i, last_par, chr) {
  # Check if it's not the sum picture but just per sample
  if (last_par != 'sum'){
    # Read file
    snps<-read.table(path,sep="\t",header=F)
    colnames(snps)<-c("gene_id", "seqnames","Start","End")
    print(head(snps))
    
    #filter_chr <- snps[snps$seqnames == chr, ]
    #print(tail(snps))
    #print(tail(filter_chr))
    
    # Join start and end separated by a '-'
    snps$range <- paste(snps$Start, "-", snps$End)
    # Make it into a GRanges, so you can make a plot with it
    genes <- GRanges(snps) #toGRanges
    # Add patient ID as a column
    snps$sampleName <- basename
    # Check whether it is the first sample or whether other samples have already been read
    if (i == 1){
      new_df <- snps
    } else {
      new_df <- rbind(df,snps)
    }
    # Return a list of genes and new_df
    return(list(genes,new_df))
  } else {
    print(head(df))
    print(tail(df))
    print('-------------')
    #df = subset(df, select=-c(sampleName))
    # Replace the gene_id column
    df$gene_id <- seq(1, nrow(df), by=1)
    # Make the columns numeric
    df[, c(1,3:4)] <- sapply(df[, c(1,3:4)], as.numeric)
    print(sum(is.na(df$Start)))
    print(sum(is.na(df$End)))
    # https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
    # Delete any NA rows
    df <- na.omit(df)
    #new_df = new_df[!is.na(new_df)]
    # Make it into a GRanges, so you can make a plot with it
    genes <- GRanges(df)
    # Return genes
    return(genes)
  }
  #, idcol = "sampleName") #, fill = TRUE) idcol = "year"
  
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
type_files <- 'other' # other or vcf
chrom <- c(1,2)


filenames <- list.files(paste("D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/", type_files, '/', sep=''), pattern="*.bed", full.names=TRUE)
print(filenames)
for (chr in chrom){
  chr <- paste('chr', chr, sep="")
  #df <- data.frame(matrix(ncol = 4, nrow = 0))
  #colnames(df) <- c("gene_id", "seqnames","Start","End")
  #colClasses = c("Date", "character", "character")
  col.names = c("gene_id", "seqnames","Start","End")
  
  df <- read.table(text = "",
                   #colClasses = colClasses,
                   col.names = col.names)
  for (i in 1:length(filenames)){
    r0=0
    r1=0.10
    if (type_files == 'other'){
      basename <- strsplit(basename(filenames[i]), "_", fixed = TRUE)[[1]][1]
    } else {
      basename <- strsplit(basename(filenames[i]), ".", fixed = TRUE)[[1]][1]
    }
    
    if ( (i%%12) == 0 || i == 1){
      if (i != 1){
        dev.off()
      }
      counter <- 0     
      
      name <- paste('D:/Hanze_Groningen/STAGE/bed/PLOTS/kary/', type_files, '/', chr, '/plot', i, basename, '.png', sep="")
      png(name, width = 2000, height = 1000) #width = 2000, height = 1000
      kp <- plotKaryotype()
      # select chromomes
      kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes=c(chr)) #chromosomes=c("chr22") main="Gene Density", 
      # add base numbers
      kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
      df_list <- create_data_karyoploteR(filenames[i], df, i, '', chr)
      
      sample1 <- df_list[[1]]
      df <- df_list[[2]]
      add <- (counter*0.14)
      kp <- add_density(kp, sample1, r0+add, r1+add, basename)
      
    } else{
      df_list <- create_data_karyoploteR(filenames[i], df, i, '', chr)
      sample1 <- df_list[[1]]
      df <- df_list[[2]]
      counter <- counter + 1
      add <- (counter*0.14)
      kp <- add_density(kp, sample1, r0+add, r1+add, basename)
      
      if (i == length(filenames)) {
        dev.off()
        name <- paste('D:/Hanze_Groningen/STAGE/bed/PLOTS/kary/', type_files, '/', chr, '/ALLplot', i, basename, '.png', sep="")
        png(name, width = 2000, height = 1000) #width = 2000, height = 1000
        kp <- plotKaryotype()
        # select chromomes
        kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes=c(chr)) #chromosomes=c("chr22") main="Gene Density", 
        # add base numbers
        kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
        sample1 <- create_data_karyoploteR(filenames[i], df, i, 'sum', chr)
        #sample1 <- df_list[[1]]
        #df <- df_list[[2]]
        # add <- (counter*0.14)
        kp <- add_density(kp, sample1, 0.0, 1.0, 'ALL')
        dev.off()        
      }
      
    }
    
  }
}


