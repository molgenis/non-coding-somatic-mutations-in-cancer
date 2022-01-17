#!/usr/bin/Rscript
library(karyoploteR)
library(GenomicRanges)

library(RSQLite)
# source("karyoploteR.R")


create_data_karyoploteR <- function(path, df, i, last_par, chr, basename_file) {
    print('---10')
    # Check if it's not the sum picture but just per sample
    if (last_par != 'sum'){
        print('---11')
        print(head(path))
        # Read file
        snps<-path
        colnames(snps)<-c("gene_id", "seqnames","Start","End")
        print(head(snps))
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
        print('---12')
        # Replace the gene_id column
        df$gene_id <- seq(1, nrow(df), by=1)
        # Make the columns numeric
        df[, c(1,3:4)] <- sapply(df[, c(1,3:4)], as.numeric)
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
    kp <- kpPlotDensity(kp, data=data, r0=r0, r1=r1, window.size = 100000)
    kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=r0, r1=r1, cex=0.8)
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
    }
    close(printer)
}


make_plots <- function(chr, filenames, df, i, basename_file, counter, r0, r1, last_par, new_add, kp) {
    print('---4')
    print(filenames)
    #
    if (new_add == 'new') {
        print('---5')
        kp <- plotKaryotype()
        # select chromomes
        kp <- plotKaryotype(genome="hg19", plot.type=1, chromosomes=c(chr)) 
        # add base numbers
        kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", cex=0.8)
        add <- (counter*0.14)
    }else{
        print('---6')
        counter <- counter + 1
        add <- (counter*0.14)
    }
    #    
    if (last_par != 'sum'){
        print('---7')
        df_list <- create_data_karyoploteR(filenames, df, i, last_par, chr, basename_file)
        print('---9')
        df <- df_list[[2]]
        
        kp <- add_density(kp, df_list[[1]], r0+add, r1+add, basename_file)
        return(list(kp, df, counter))
    }else{
        print('---8')
        all_snps <- create_data_karyoploteR(filenames, df, i, last_par, chr, basename_file)
        kp <- add_density(kp, all_snps, 0.0, 1.0, 'ALL')
        return(kp)
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

    # Parameters for where the plots will be placed
    r0=0
    r1=0.10

    per_project = FALSE #TRUE or FALSE
    

    conn <- dbConnect(RSQLite::SQLite(), db_path)
    # num_projects <- dbGetQuery(conn, "SELECT ID, chr, pos_start, pos_end FROM snp")
    projects <- dbGetQuery(conn, "SELECT ID, project_ID  FROM project")

    # for (chr in chrom){
    i = 0
    # Make empty dataframe
    df <- read.table(text = "", col.names = c("gene_id", "seqnames","Start","End"))
    for (project_ID in nrow(projects)){ #
        ID_project <- projects['ID'][[1]][project_ID]
        basename_project <- projects['project_ID'][[1]][project_ID]
        if (!dir.exists(path_info_save)){
            dir.create(path_info_save, recursive = T)
        } 
        


        print(project_ID)
        donors <- dbGetQuery(conn, sprintf("SELECT ID, donor_ID FROM donor WHERE project_ID =%s", ID_project))
        for (donor_ID in nrow(donors)){ #['ID'][[1]]
            i = i + 1
            print(donor_ID)
            ID_donor <- donors['ID'][[1]][donor_ID]
            basename_donor <- donors['donor_ID'][[1]][donor_ID]
            donor_snp_id <- dbGetQuery(conn, sprintf("SELECT snp_ID FROM donor_has_snp WHERE donor_project_ID =%s AND donor_ID = %s", ID_project, ID_donor))
            #snp_df = data.frame(matrix("", ncol = 4, nrow = nrow(donor_snp_id)))
            #colnames(snp_df)<-c("ID", "chr","pos_start","pos_end")    
            snp_df <- read.table(text = "", col.names = c("ID", "chr","pos_start","pos_end"))
            #for (snp_ID in donor_snp_id['snp_ID'][[1]]){
            print('TEST1')
            for (snp_ID in seq(1, nrow(donor_snp_id), by=1)){
                ID_snp = donor_snp_id['snp_ID'][[1]][snp_ID]
                print(ID_snp)
                # print(basename_donor)
                # print(dbGetQuery(conn, "SELECT count(ID) FROM snp WHERE ID=%s AND seq_strategy = 'WGS' AND chr = %s"))
                the_snps <- dbGetQuery(conn, sprintf("SELECT ID, chr, pos_start, pos_end FROM snp WHERE ID=%s AND seq_strategy = 'WGS'", ID_snp)) #AND chr = %s
                # print('nooo')
                snp_df <- rbind(snp_df, the_snps)
                #snp_df[snp_ID,] = test
                }
            }
            for (chr in chrom){
                #TODO
                # TODO
            }
            print('TEST2')
            print(nrow(snp_df))
            if (nrow(snp_df) > 0){
                print(head(snp_df))
                # paste 'chr' for each chromosome
                chr <- paste('chr', chr, sep="")
                print('hoi')
                # Checked that a maximum of 12 participants are clearly visible in a plot. 
                # So now he makes multiple plots with 12 (or fewer) participants per plot.
                print('----------')
                if ( (i%%12) == 0 || i == 1){
                    print('---1')
                    if (i != 1){
                    dev.off()
                    }
                    counter <- 0            
                    name <- paste(path_info_save, 'plot', i, '.png', sep="")
                    png(name, width = 2000, height = 1000) 
                    info <- make_plots(chr, snp_df, df, i, basename_donor, counter, r0, r1, '', 'new', '') #kp, df, counter
                    kp <- info[[1]]
                    df <- info[[2]]
                    counter <- info[[3]]        
                } else{
                    print('---2')
                    info <- make_plots(chr, snp_df, df, i, basename_donor, counter, r0, r1, '', 'add', kp) #kp, df, counter
                    kp <- info[[1]]
                    df <- info[[2]]
                    counter <- info[[3]]  
                    # When you arrive at the last *.bed file, a plot is made in which all patients are added up. 
                    # So that you have everything in 1 enumerating picture.
                    if (i == nrow(snp_df)) {
                        print('---3')
                        dev.off()
                        name <- paste(path_info_save, 'ALLplot_chr', chr, '.png', sep="")
                        peaks_file <- paste(path_info_save, 'ALLplot_chr', chr, '.tsv', sep="")
                        png(name, width = 2000, height = 1000)
                        kp <- make_plots(chr, snp_df, df, i, basename_donor, counter, r0, r1, 'sum', 'new', '')
                        peaks(kp, peaks_file)
                        dev.off()        
                    }        
                }    

            }else{
                print('LEEG')
            }            
    }
    # }    
}

main()



