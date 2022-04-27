#!/usr/bin/Rscript
library(karyoploteR)
library(GenomicRanges)
library(RSQLite)
#library(zeallot) #https://stackoverflow.com/questions/13538628/returning-different-data-frames-in-a-function-r



create_data_karyoploteR <- function(path, df, i, last_par, chr, basename_file) {
    ##########################################
    #
    ##########################################
    # Check if it's not the sum picture but just per sample
    if (last_par != 'sum'){
          # Read file
          snps<-read.table(path,sep="\t",header=F)
          colnames(snps)<-c("gene_id", "seqnames","Start","End")
          snpsHOI <- snps[snps$seqnames == paste("chr", chr),]
          if (dim(snpsHOI)[1] == 0) {
              print('JAAAAA')
          }else{
              print('neee')
          }
          # Join start and end separated by a '-'
          snps$range <- paste(snps$Start, "-", snps$End)
          # Make it into a GRanges, so you can make a plot with it
          all_snps <- GRanges(snps) #toGRanges
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
          print(head(df))
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
    #
    ##########################################
    # Plot the SNP density over the genome
    kp <- kpPlotDensity(kp, data=data, r0=r0, r1=r1, window.size = 100000)
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
    #
    ##########################################
    # Get information out of the plots
    density <- kp$latest.plot$computed.values$density
    windows <- kp$latest.plot$computed.values$windows
    # Make file for the information
    textfile=file.path(peaks_file);    
    printer = file(textfile,"a+");
    # Get the 10 highest peaks
    indexes <- order(density, decreasing=TRUE)[1:20]
    count <- 0
    # Grab the window belonging to the highest peaks
    for (index in indexes){
        count <- count + 1
        # Write info to file - 'hight of peak' - 'start of window' - 'end of window'
        #write(c(density[index], start(windows[index]), end(windows[index])), textfile,sep = "\t",append = TRUE, ncolumns = 3);
        write(c(count, chr, start(windows[index]), end(windows[index])), textfile,sep = "\t",append = TRUE, ncolumns = 4);
    }
    close(printer)
}


make_plots <- function(chr, filenames, df, i, basename_file, counter, r0, r1, last_par, new_add, kp, num_for_r, custom.cytobands, zoom.region) {
    ##########################################
    #
    ##########################################
    # With a completely new plot, an object has to be created first. 
    # For a plot that has already been created, only a plot needs to be added.
    if (new_add == 'new') {
        # Create karyoplot object
        #kp <- plotKaryotype()
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
        df_list <- create_data_karyoploteR(filenames[i], df, i, last_par, chr, basename_file)
        # The data frame that is continuously expanded with new data for the SUM plot
        df <- df_list[[2]]
        # Call add_density
        kp <- add_density(kp, df_list[[1]], r0+add, r1+add, basename_file)
        return(list(kp, df, counter))
    }else{
        # Call create_data_karyoploteR
        all_snps <- create_data_karyoploteR(filenames[i], df, i, last_par, chr, basename_file)
        # Call add_density
        kp <- add_density(kp, all_snps, 0.0, 1.0, 'ALL')
        return(kp)
    }
}


read_files <- function(filenames, chr, r0, r1, path_info_save, path_peaks_chr, num_of_pictures, num_for_r, custom.cytobands) {
    ##########################################
    #
    ##########################################
    print(filenames)
    # Paste 'chr' for each chromosome
    chr <- paste('chr', chr, sep="")
    # Make empty dataframe
    df <- read.table(text = "", col.names = c("gene_id", "seqnames","Start","End"))
    # Loop over all files that end with .bed
    for (i in 1:length(filenames)){  
        print(i)
        # With one file I have to split on something else than with the other file
        basename_file <- tools::file_path_sans_ext(strsplit(basename(filenames[i]), ".", fixed = TRUE)[[1]][1])
        peaks_file_big <- paste(path_peaks_chr, basename_file, '_', chr, '.tsv', sep="")
        print(peaks_file_big)
        peaks <- read.table(file = peaks_file_big, sep = '\t', header = FALSE, col.names = c("id", "chromosoom","Start","End"))
        for (index_row in 1:nrow(peaks)){
            row <- peaks[index_row,]
            zoom.region <- toGRanges(data.frame(row$chromosoom, row$Start, row$End))
            print(basename_file)
            # Checked that a maximum of 12 participants are clearly visible in a plot. 
            # So now he makes multiple plots with 12 (or fewer) participants per plot.
            if ( (i%%(num_of_pictures+1)) == 0 || i == 1){
                # If i is not 1, he has to close the previous plot, with 1 that is not yet possible,
                # because then no plot has been made yet.
                if (i != 1 || num_of_pictures == 1){
                    dev.off()
                }
                # Set counter to 0: Keeps track of plot number
                counter <- 0 
                # Create name for the figure that saves ~12 plots      
                name <- paste(path_info_save, 'plot', i, '_', basename_file, '_', row$id, '_', row$chromosoom, '_', row$Start, '_', row$End, '.png', sep="")
                # Make figure
                png(name, width = 2000, height = 1000) 
                # Call make_plots
                info <- make_plots(chr, filenames, df, i, basename_file, counter, r0, r1, '', 'new', '', num_for_r, custom.cytobands, zoom.region) #kp, df, counter
                # The kp plot (karyoplot object)
                kp <- info[[1]]
                # The data frame that is continuously expanded with new data for the SUM plot
                df <- info[[2]]
                # counter: Keeps track of plot number
                counter <- info[[3]]  
                peaks_file <- paste(path_info_save, basename_file, '_', chr, '_', row$id, '_', row$chromosoom, '_', row$Start, '_', row$End, '.tsv', sep="")
                peaks(kp, peaks_file, chr)
            } else{
                # Call make_plots
                info <- make_plots(chr, filenames, df, i, basename_file, counter, r0, r1, '', 'add', kp, num_for_r, custom.cytobands, zoom.region) #kp, df, counter
                # The kp plot (karyoplot object)
                kp <- info[[1]]
                # The data frame that is continuously expanded with new data for the SUM plot
                df <- info[[2]]
                # counter: Keeps track of plot number
                counter <- info[[3]]  
                peaks_file <- paste(path_info_save, basename_file, '_', chr, '_', row$id, '_', row$chromosoom, '_', row$Start, '_', row$End, '.tsv', sep="")
                peaks(kp, peaks_file, chr)
                # When you arrive at the last *.bed file, a plot is made in which all patients are added up. 
                # So that you have everything in 1 enumerating picture.
                if (i == length(filenames)){
                    dev.off()
                    # Create name for the figure that saves SUM of plots
                    name <- paste(path_info_save, 'ALLplot_', i, '_', basename_file, '_', row$id, '_', row$chromosoom, '_', row$Start, '_', row$End, '.png', sep="")
                    peaks_file <- paste(path_info_save, 'ALLplot_', chr, '_', row$id, '_', row$chromosoom, '_', row$Start, '_', row$End, '.tsv', sep="")
                    png(name, width = 2000, height = 1000)
                    # Call make_plots
                    kp <- make_plots(chr, filenames, df, i, basename_file, counter, r0, r1, 'sum', 'new', '', num_for_r, custom.cytobands, zoom.region)
                    # Call peaks
                    peaks(kp, peaks_file, chr)
                    dev.off()     
                }   
            }
        }       
    }   
}









read_files <- function(filenames, chr, r0, r1, path_info_save, path_peaks_chr, num_of_pictures, num_for_r, custom.cytobands, zoom.region, row) {
    ##########################################
    #
    ##########################################
    print(filenames)
    # Paste 'chr' for each chromosome
    chr <- paste('chr', chr, sep="")
    # Make empty dataframe
    df <- read.table(text = "", col.names = c("gene_id", "seqnames","Start","End"))
    # Loop over all files that end with .bed
    for (i in 1:length(filenames)){  
        print(i)
        # With one file I have to split on something else than with the other file
        basename_file <- tools::file_path_sans_ext(strsplit(basename(filenames[i]), ".", fixed = TRUE)[[1]][1])
        # peaks_file_big <- paste(path_peaks_chr, basename_file, '_', chr, '.tsv', sep="")
        # print(peaks_file_big)
        # peaks <- read.table(file = peaks_file_big, sep = '\t', header = FALSE, col.names = c("id", "chromosoom","Start","End"))
        # for (index_row in 1:nrow(peaks)){
        # row <- peaks[index_row,]
        # zoom.region <- toGRanges(data.frame(row$chromosoom, row$Start, row$End))
        print(basename_file)
        # Checked that a maximum of 12 participants are clearly visible in a plot. 
        # So now he makes multiple plots with 12 (or fewer) participants per plot.
        if ( (i%%(num_of_pictures+1)) == 0 || i == 1){
            # If i is not 1, he has to close the previous plot, with 1 that is not yet possible,
            # because then no plot has been made yet.
            if (i != 1 || num_of_pictures == 1){
                dev.off()
            }
            # Set counter to 0: Keeps track of plot number
            counter <- 0 
            # Create name for the figure that saves ~12 plots      
            name <- paste(path_info_save, 'plot', i, '_', basename_file, '_', row$id, '_', row$chromosoom, '_', row$Start, '_', row$End, '.png', sep="")
            # Make figure
            png(name, width = 2000, height = 1000) 
            # Call make_plots
            info <- make_plots(chr, filenames, df, i, basename_file, counter, r0, r1, '', 'new', '', num_for_r, custom.cytobands, zoom.region) #kp, df, counter
            # The kp plot (karyoplot object)
            kp <- info[[1]]
            # The data frame that is continuously expanded with new data for the SUM plot
            df <- info[[2]]
            # counter: Keeps track of plot number
            counter <- info[[3]]  
            peaks_file <- paste(path_info_save, basename_file, '_', chr, '_', row$id, '_', row$chromosoom, '_', row$Start, '_', row$End, '.tsv', sep="")
            peaks(kp, peaks_file, chr)
        } else{
            # Call make_plots
            info <- make_plots(chr, filenames, df, i, basename_file, counter, r0, r1, '', 'add', kp, num_for_r, custom.cytobands, zoom.region) #kp, df, counter
            # The kp plot (karyoplot object)
            kp <- info[[1]]
            # The data frame that is continuously expanded with new data for the SUM plot
            df <- info[[2]]
            # counter: Keeps track of plot number
            counter <- info[[3]]  
            peaks_file <- paste(path_info_save, basename_file, '_', chr, '_', row$id, '_', row$chromosoom, '_', row$Start, '_', row$End, '.tsv', sep="")
            peaks(kp, peaks_file, chr)
            # When you arrive at the last *.bed file, a plot is made in which all patients are added up. 
            # So that you have everything in 1 enumerating picture.
            if (i == length(filenames)){
                dev.off()
                # Create name for the figure that saves SUM of plots
                name <- paste(path_info_save, 'ALLplot_', i, '_', basename_file, '_', row$id, '_', row$chromosoom, '_', row$Start, '_', row$End, '.png', sep="")
                peaks_file <- paste(path_info_save, 'ALLplot_', chr, '_', row$id, '_', row$chromosoom, '_', row$Start, '_', row$End, '.tsv', sep="")
                png(name, width = 2000, height = 1000)
                # Call make_plots
                kp <- make_plots(chr, filenames, df, i, basename_file, counter, r0, r1, 'sum', 'new', '', num_for_r, custom.cytobands, zoom.region)
                # Call peaks
                peaks(kp, peaks_file, chr)
                dev.off()     
            }   
        }
        # }       
    }   
}



main <- function() {
    # List of chromosomes
    #numbers <-  seq(1, 22, by=1)
    #chrom <- append(numbers, c('X', 'Y'))
    chrom <- c(3)
    
    
    # Cancer
    #path_file = "D:/Hanze_Groningen/STAGE/R/select_inf/cancer/"
    #path_save = "D:/Hanze_Groningen/STAGE/R/PLOTS/kary/cancer/"
    # Project
    #path_file = "D:/Hanze_Groningen/STAGE/R/select_inf/project/"
    #path_save = "D:/Hanze_Groningen/STAGE/R/PLOTS/kary/project/"
    # VS
    path_file = "D:/Hanze_Groningen/STAGE/R/select_inf/vs/"
    path_save = "D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/ZOOM/"
    
    path_peaks = "D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/before/"
    
    
    # Color chromosome on genes
    custom.cytobands <- toGRanges('D:/Hanze_Groningen/STAGE/db/mycytobands_R.txt')
    # How many layers do you want on one plot
    num_of_pictures = 2
    # Number for calculating R
    num_for_r = 12/num_of_pictures
    # Parameters for where the plots will be placed
    r0=0 * num_for_r
    r1=0.10 * num_for_r
    
    
    # If per_project is TRUE, you need to connect to the database to get the different project names.
    # Get all the filenames in that path that ends with .bed
    filenames <- list.files(path_file, pattern="*.tsv", full.names=TRUE)
    # Delete files with file.size > 1
    filenames <- filenames[sapply(filenames, file.size) > 1]
    
    # Loop over chromosomes
    for (chr in chrom){
        # peaks_file_big <- paste(path_peaks_chr, basename_file, '_', chr, '.tsv', sep="")
        # print(peaks_file_big)
        peaks <- read.table(file = 'D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/before/chr3/breast_chr3.tsv', sep = '\t', header = FALSE, col.names = c("id", "chromosoom","Start","End"))
        for (index_row in 1:nrow(peaks)){
            row <- peaks[index_row,]
            zoom.region <- toGRanges(data.frame(row$chromosoom, row$Start, row$End))
            print(chr)
            path_info_save = paste(path_save, '/chr', chr, '/', sep='')
            path_peaks_chr = paste(path_peaks, '/chr', chr, '/', sep='')
            print(path_info_save)
            # Makes folder if it doesn't exists
            if (!dir.exists(path_info_save)){
                dir.create(path_info_save, recursive = T)
            }      
            # Call read_files      
            read_files(filenames, chr, r0, r1, path_info_save, path_peaks_chr, num_of_pictures, num_for_r, custom.cytobands, zoom.region, row)  
        }      
    }
}

# runs only when script is run by itself
if (sys.nframe() == 0){
    main()
}