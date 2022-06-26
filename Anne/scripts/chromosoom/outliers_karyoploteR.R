#!/usr/bin/Rscript
library(karyoploteR)
library(GenomicRanges)



create_data_karyoploteR <- function(snps, df, i, last_par, chr, basename_file) {
    ##########################################
    # Create karyoplote data for the plot
    ##########################################
    # Check if it's not the sum picture but just per sample
    if (last_par != 'sum'){
        # Join start and end separated by a '-'
        snps$range <- paste(snps$Start, "-", snps$End)
        # Make it into a GRanges, so you can make a plot with it
        all_snps <- GRanges(snps) 
        # Add patient ID as a column
        snps$sampleName <- basename_file
        # Check whether it is the first sample or whether other samples have already been read
        if (i == 1){
            new_df <- snps
        } else {
            new_df <- rbind(df,snps)
        }
        # Return a list of all_snps and new_df
        return(list(all_snps,new_df))
    } else {
        # Replace the gene_id column
        df$gene_id <- seq(1, nrow(df), by=1)
        # Make the columns numeric
        df[, c(1,3:4)] <- sapply(df[, c(1,3:4)], as.numeric)
        # Delete any NA rows
        df <- na.omit(df)
        # Make it into a GRanges, so you can make a plot with it
        all_snps <- GRanges(df)
        # Return all_snps
        return(all_snps)
    }  
}

add_density <- function(kp, data, r0, r1, label) {
    ##########################################
    # Add density, axis and labels to the plot
    ##########################################
    # Plot the SNP density over the genome
    kp <- kpPlotDensity(kp, data=data, r0=r0, r1=r1, window.size = 100)
    # Adding Y axis to our plots may help interpreting the magnitudes of plotted data and
    # identifying the space dedicated to each plot.
    kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1, cex=0.8)
    # Add labels. Labels will appear on the left side of the plot.
    kpAddLabels(kp, labels=label, label.margin=0.04, r0=r0, r1=r1, data.panel = 1, cex=1.2, col = "cornflowerblue")
    # Return karyoplot object
    return(kp)
}

peaks <- function(kp, peaks_file, chr){
    ##########################################
    # Get the peaks
    ##########################################
    # Get information out of the plots
    density <- kp$latest.plot$computed.values$density
    windows <- kp$latest.plot$computed.values$windows
    # Make file for the information
    textfile=file.path(peaks_file);    
    printer = file(textfile,"a+");
    # Get the 10 highest peaks
    indexes <- order(density, decreasing=TRUE)#[1:20]
    count <- 0
    # Grab the window belonging to the highest peaks
    for (index in indexes){
        count <- count + 1
        # Write info to file - 'hight of peak' - 'start of window' - 'end of window'
        write(c(count, density[index], chr, start(windows[index]), end(windows[index])), textfile,sep = "\t",append = TRUE, ncolumns = 5);
    }
    close(printer)
}


make_plots <- function(chr, filenames, df, i, basename_file, counter, r0, r1, last_par, new_add, kp, num_for_r, custom.cytobands, zoom.region) {
    ##########################################
    # Create the plots
    ##########################################
    # With a completely new plot, an object has to be created first. 
    # For a plot that has already been created, only a plot needs to be added.
    if (new_add == 'new') {
        # Select chromosome
        kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes=c(chr), cytobands = custom.cytobands, zoom=zoom.region) 
        # Add base numbers: add the base numbering below each chromosome
        kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 10000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
        # add: How much should be added in terms of margin for the next plot
        add <- (counter*0.14*num_for_r)
    }else{
        # counter: Keeps track of plot number
        counter <- counter + 1
        # add: How much should be added in terms of margin for the next plot
        add <- (counter*0.14*num_for_r)
    }
    # For all plots above each other (entire plot) other things have to happen than with the SUM of all plots
    if (last_par != 'sum'){
        # Call create_data_karyoploteR
        df_list <- create_data_karyoploteR(filenames, df, i, last_par, chr, basename_file)
        # The data frame that is continuously expanded with new data for the SUM plot
        df <- df_list[[2]]
        # Call add_density
        kp <- add_density(kp, df_list[[1]], r0+add, r1+add, basename_file)
        return(list(kp, df, counter))
    }else{
        # Call create_data_karyoploteR
        all_snps <- create_data_karyoploteR(filenames, df, i, last_par, chr, basename_file)
        # Call add_density
        kp <- add_density(kp, all_snps, 0.0, 1.0, 'ALL')
        return(kp)
    }
}


set_plot_info <- function(select_breast, select_nonbreast, chr, start, end, r0, r1, path_info_save, num_of_pictures, num_for_r, custom.cytobands, zoom.region) {
    ##########################################
    # 
    ##########################################
    i <- 1
    print(i)
    # Checked that a maximum of 12 participants are clearly visible in a plot. 
    # So now he makes multiple plots with 12 (or fewer) participants per plot.
    basename_file <- 'breast'
    # Set counter to 0: Keeps track of plot number
    counter <- 0 
    # Create name for the figure that saves ~12 plots      
    name <- paste(path_info_save, 'plot', chr,'_start_', start,'_end_', end, '.png', sep="")
    # Make figure
    png(name, width = 2000, height = 1000) 
    # Call make_plots
    info <- make_plots(chr, select_breast, df, i, basename_file, counter, r0, r1, '', 'new', '', num_for_r, custom.cytobands, zoom.region) #kp, df, counter
    # The kp plot (karyoplot object)
    kp <- info[[1]]
    # The data frame that is continuously expanded with new data for the SUM plot
    df <- info[[2]]
    # counter: Keeps track of plot number
    counter <- info[[3]]  
    dev.off()
        
}




path_outliers = "D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/before/outliers_unique_breast.tsv"
path_breast = "D:/Hanze_Groningen/STAGE/R/select_inf/vs/breast.tsv"
path_nonbreast = "D:/Hanze_Groningen/STAGE/R/select_inf/vs/nonbreast.tsv"
outliers <- read.table(file = path_outliers, sep = '\t', header = TRUE)[1,]
breast <- read.table(file = path_breast, sep = '\t', header = FALSE)
colnames(breast) <- c("gene_id", "seqnames","Start","End")
nonbreast <- read.table(file = path_nonbreast, sep = '\t', header = FALSE)
colnames(nonbreast) <- c("gene_id", "seqnames","Start","End")
path_info_save = 'D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/ZOOM/outliers/'

# Color chromosome on genes
custom.cytobands <- toGRanges('D:/Hanze_Groningen/STAGE/db/mycytobands_R.txt')
# How many layers do you want on one plot
num_of_pictures = 2
# Number for calculating R
num_for_r = 12/num_of_pictures
# Parameters for where the plots will be placed
r0=0 * num_for_r
r1=0.10 * num_for_r

for (i in 1:nrow(outliers)){
  outl <- outliers[i,]
  print(outl)
  select_breast <- breast[breast$Start>=outl$start_region & breast$Start<=outl$end_region & breast$seqnames==outl$chr,]
  select_nonbreast <- nonbreast[nonbreast$Start>=(outl$start_region+1000000) & nonbreast$Start<=(outl$end_region+1000000) & nonbreast$seqnames==outl$chr, ]
  zoom.region <- toGRanges(data.frame(outl$chr, (outl$start_region+1000000), (outl$end_region+1000000)))
  set_plot_info(select_breast, select_nonbreast, outl$chr, (outl$start_region+1000000), (outl$end_region+1000000), r0, r1, path_info_save, num_of_pictures, num_for_r, custom.cytobands, zoom.region)
  
}



