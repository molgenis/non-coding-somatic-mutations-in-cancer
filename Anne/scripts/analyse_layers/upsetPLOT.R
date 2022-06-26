library(UpSetR)
vector_all_afterGene <- read.table(file = "E:/analyse/new/correction/snps/all_list/ALL_afterGene.tsv", sep = '\t', header = TRUE)[,1]
vector_all_beforeGene <- read.table(file = "E:/analyse/new/correction/snps/all_list/ALL_beforeGene.tsv", sep = '\t', header = TRUE)[,1]
vector_all_DNase <- read.table(file = "E:/analyse/new/correction/snps/all_list/ALL_DNase.tsv", sep = '\t', header = TRUE)[,1]
vector_all_persnp <- read.table(file = "E:/analyse/new/correction/snps/all_list/ALL_persnp.tsv", sep = '\t', header = TRUE)[,1]
vector_all_TFBS <- read.table(file = "E:/analyse/new/correction/snps/all_list/ALL_TFBS.tsv", sep = '\t', header = TRUE)[,1]
vector_all_UCNE <- read.table(file = "E:/analyse/new/correction/snps/all_list/ALL_UCNE.tsv", sep = '\t', header = TRUE)[,1]

sigs_list_all <- list('afterGene' = vector_all_afterGene, 'beforeGene' = vector_all_beforeGene, 'DNase'=vector_all_DNase, 
                  'persnp' = vector_all_persnp, 'TFBS' = vector_all_TFBS, 'UCNE'=vector_all_UCNE)

upset(fromList(sigs_list_all), order.by = 'freq', nsets = 6) #queries = queries, 