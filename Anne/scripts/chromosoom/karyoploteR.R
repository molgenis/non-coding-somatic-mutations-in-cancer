#!/usr/bin/Rscript
library(karyoploteR)
library(GenomicRanges)
library(RSQLite)
#library(zeallot) #https://stackoverflow.com/questions/13538628/returning-different-data-frames-in-a-function-r



create_data_karyoploteR <- function(path, df, i, last_par, chr, basename_file) {
    # Check if it's not the sum picture but just per sample
    if (last_par != 'sum'){
        # Read file
        snps<-read.table(path,sep="\t",header=F)
        colnames(snps)<-c("gene_id", "seqnames","Start","End")    
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
        # Replace the gene_id column
        df$gene_id <- seq(1, nrow(df), by=1)
        # Make the columns numeric
        df[, c(1,3:4)] <- sapply(df[, c(1,3:4)], as.numeric)
        # print(sum(is.na(df$Start)))
        # print(sum(is.na(df$End)))
        # https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
        # Delete any NA rows
        df <- na.omit(df)
        # Make it into a GRanges, so you can make a plot with it
        all_snps <- GRanges(df)
        # Return all_snps
        return(all_snps)
    }  
}

add_density <- function(kp, data, r0, r1, label) {
    # plot the SNP density over the genome
    kp <- kpPlotDensity(kp, data=data, r0=r0, r1=r1, window.size = 100000) #window.size = 0.5e6 100000
    kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1, cex=0.8)
    #kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1)
    kpAddLabels(kp, labels=label, label.margin=0.04, r0=r0, r1=r1, data.panel = 1, cex=1.2, col = "cornflowerblue")
    return(kp)
}

peaks <- function(kp, peaks_file){
    density <- kp$latest.plot$computed.values$density
    windows <- kp$latest.plot$computed.values$windows
    textfile=file.path(peaks_file);
    printer = file(textfile,"a+");
    indexes <- order(density, decreasing=TRUE)[1:10]
    for (index in indexes){
        # https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GRanges_and_GRangesList_slides.pdf
        write(c(density[index], start(windows[index]), end(windows[index])), textfile,sep = "\t",append = TRUE, ncolumns = 3);
        #write("\n", textfile, append=TRUE)
    }
    close(printer)
}


make_plots <- function(chr, filenames, df, i, basename_file, counter, r0, r1, last_par, new_add, kp) {
    #
    if (new_add == 'new') {
        kp <- plotKaryotype()
        # select chromomes
        kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes=c(chr)) 
        # add base numbers
        kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
        add <- (counter*0.14)
    }else{
        counter <- counter + 1
        add <- (counter*0.14)
    }
    #    
    if (last_par != 'sum'){
        df_list <- create_data_karyoploteR(filenames[i], df, i, last_par, chr, basename_file)
        df <- df_list[[2]]
        
        kp <- add_density(kp, df_list[[1]], r0+add, r1+add, basename_file)
        return(list(kp, df, counter))
    }else{
        all_snps <- create_data_karyoploteR(filenames[i], df, i, last_par, chr, basename_file)
        kp <- add_density(kp, all_snps, 0.0, 1.0, 'ALL')
        return(kp)
    }
}


read_files <- function(type_files, filenames, chr, r0, r1, path_info_save) {
    # paste 'chr' for each chromosome
    chr <- paste('chr', chr, sep="")
    # Make empty dataframe
    df <- read.table(text = "", col.names = c("gene_id", "seqnames","Start","End"))
    # Loop over all files that end with .bed
    for (i in 1:length(filenames)){        
        # With one file I have to split on something else than with the other file
        if (type_files == 'other'){
            basename_file <- tools::file_path_sans_ext(strsplit(basename(filenames[i]), "_", fixed = TRUE)[[1]][2])
        } else {
            basename_file <- strsplit(basename(filenames[i]), ".", fixed = TRUE)[[1]][1]
        }
        # Checked that a maximum of 12 participants are clearly visible in a plot. 
        # So now he makes multiple plots with 12 (or fewer) participants per plot.
        if ( (i%%12) == 0 || i == 1){
            if (i != 1){
            dev.off()
            }
            counter <- 0            
            name <- paste(path_info_save, 'plot', i, '.png', sep="")
            png(name, width = 2000, height = 1000) 
            info <- make_plots(chr, filenames, df, i, basename_file, counter, r0, r1, '', 'new', '') #kp, df, counter
            kp <- info[[1]]
            df <- info[[2]]
            counter <- info[[3]]        
        } else{
            info <- make_plots(chr, filenames, df, i, basename_file, counter, r0, r1, '', 'add', kp) #kp, df, counter
            kp <- info[[1]]
            df <- info[[2]]
            counter <- info[[3]]  
            # When you arrive at the last *.bed file, a plot is made in which all patients are added up. 
            # So that you have everything in 1 enumerating picture.
            if (i == length(filenames)) {
                dev.off()
                name <- paste(path_info_save, 'ALLplot_chr', chr, '.png', sep="")
                peaks_file <- paste(path_info_save, 'ALLplot_chr', chr, '.tsv', sep="")
                png(name, width = 2000, height = 1000)
                kp <- make_plots(chr, filenames, df, i, basename_file, counter, r0, r1, 'sum', 'new', '')
                peaks(kp, peaks_file)
                dev.off()        
            }        
        }    
    }
}



main <- function() {
    # List of chromosomes
    numbers <-  seq(1, 2, by=1)
    chrom <- append(numbers, c('X', 'Y'))
    # Which type of file
    type_files <- 'other' # other or vcf

    path_file = "D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/"
    path_save = "D:/Hanze_Groningen/STAGE/bed/PLOTS/kary/"
    db_path <- "D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene.db"

    per_project = FALSE #TRUE or FALSE

    # Parameters for where the plots will be placed
    r0=0
    r1=0.10

    if (per_project){
        # Create a connection to our new database, CarsDB.db
        # you can check that the .db file has been created on your working directory
        conn <- dbConnect(RSQLite::SQLite(), db_path)
        num_projects <- dbGetQuery(conn, "SELECT project_ID FROM project")
        for (project_ID in seq(1, nrow(num_projects), by=1)){
            # print(num_projects[project_ID,1])
            # Get all the filenames in that path that ends with .bed
            patt = sprintf("%s_*.bed", project_ID)
            filenames <- list.files(paste(path_file, type_files, '/', sep=''), pattern=glob2rx(patt), full.names=TRUE)
            # delete files with file.size > 1
            filenames <- filenames[sapply(filenames, file.size) > 1]
            print(length(filenames))
            for (chr in chrom){
                print(chr)
                path_info_save = paste(path_save, type_files, '/', chr, '/', num_projects[project_ID,1],'/', sep='')
                if (!dir.exists(path_info_save)){
                    dir.create(path_info_save, recursive = T)
                }  
                read_files(type_files, filenames, chr, r0, r1, path_info_save)        
            }

        }
    }else{
        # Get all the filenames in that path that ends with .bed
        filenames <- list.files(paste(path_file, type_files, '/', sep=''), pattern="*.bed", full.names=TRUE)  
        # delete files with file.size > 1
        filenames <- filenames[sapply(filenames, file.size) > 1]   
        # Loop over chromsomes
        for (chr in chrom){
            path_info_save = paste(path_save, type_files, '/chr', chr, '/ALL/', sep='')
            print(path_info_save)
            if (!dir.exists(path_info_save)){
                dir.create(path_info_save, recursive = T)
            }            
            read_files(type_files, filenames, chr, r0, r1, path_info_save)        
        }
    }    
}



#https://stackoverflow.com/questions/2968220/is-there-an-r-equivalent-of-the-pythonic-if-name-main-main
# main()
# if (interactive()) {
#   main()
# }

# runs only when script is run by itself
if (sys.nframe() == 0){
    main()
}