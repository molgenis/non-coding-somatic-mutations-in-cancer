#!/usr/bin/env Rscript
library(chromoMap)

# chromosome files
chr_file = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/length_chromosome.tsv"
# annotation files
anno_file = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/data/SS6004099_merge_manual_bowtie.tsv"

head(read.table(chr_file,sep = "\t"))
head(read.table(anno_file,sep = "\t"))

chr = chromoMap(chr_file,anno_file, chr_width = 50, chr_length = 15,  canvas_width = 5000,
                canvas_height = 1000)
# Save an object to a file
saveRDS(chr, file = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/my_data.rds")
#save.image(file = "my_work_space.rds")

#saveWidget(chr, "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/p1.html")

saveWidget(widget=chr, file="/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/p2.html", selfcontained=TRUE, knitrOptions=list())


chr2 = chromoMap(chr_file,anno_file,
                 data_based_color_map = T,
                 data_type = "numeric",
                 plots = "scatter", chr_width = 50, 
                 chr_length = 15,  canvas_width = 5000,
                 canvas_height = 1000)

saveWidget(widget=chr2, file="/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/p22.html")#, selfcontained=TRUE, knitrOptions=list())


chr3 = chromoMap(chr_file,anno_file,
                 data_based_color_map = T,
                 data_type = "numeric",
                 plots = "tags")

saveWidget(widget=chr3, file="/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/p3.html")#, selfcontained=TRUE, knitrOptions=list())

chr4 = chromoMap(chr_file,anno_file,
                 data_based_color_map = T,
                 data_type = "numeric",
                 plots = "bar",
                 heat_map = F)
saveWidget(widget=chr4, file="/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/p4.html")#, selfcontained=TRUE, knitrOptions=list())