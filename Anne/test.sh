#!/usr/bin/bash

#samtools view -h SS6005044.sorted.bam chr22 > chr22/5044/SS6005044.chr22.sam
#cd chr22/5044
#samtools view -bS SS6005044.chr22.sam > SS6005044.chr22.bam
#samtools sort SS6005044.chr22.bam -o SS6005044_chr22.sort.bam
#samtools sort -n SS6005044.chr22.bam -o SS6005044_name_chr22.sort.bam


cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/5044

ml Anaconda3/5.3.0
ml SAMtools/1.9-foss-2018b
ml GATK/4.1.4.1-Java-8-LTS
ml BCFtools/1.11-GCCcore-7.3.0
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS


#reads that mapped properly as pairs
samtools view -u -f 1 -F 12 SS6005044_chr22.sort.bam > SS6005044_chr22_map_map.bam
# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 SS6005044_chr22.sort.bam > SS6005044_chr22_unmap_map.bam
# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 SS6005044_chr22.sort.bam > SS6005044_chr22_map_unmap.bam
# R1 & R2 unmapped
samtools view -u -f 12 -F 256 SS6005044_chr22.sort.bam > SS6005044_chr22_unmap_unmap.bam
# merge the three files
samtools merge -u SS6005044_chr22_unmapped.bam SS6005044_chr22_unmap_map.bam SS6005044_chr22_map_unmap.bam SS6005044_chr22_unmap_unmap.bam
# BAM files must be resorted so that they are ordered by read ID instead of location in the reference
samtools sort -n SS6005044_chr22_map_map.bam -o SS6005044_chr22_mapped.sort.bam
samtools sort -n SS6005044_chr22_unmapped.bam -o SS6005044_chr22_unmapped.sort.bam
# it is a good idea to check that you have the correct number of reads and no redundancy. You can summarize the original BAM file to get an idea of where you started.
samtools flagstat SS6005044_chr22.sort.bam
samtools view -c SS6005044_chr22_mapped.sort.bam
samtools view -c SS6005044_chr22_unmapped.sort.bam
# extract the FASTQ reads into two paired read files
#    X bamToFastq -i SS6005044_chr22_mapped.sort.bam -fq SS6005044_chr22_mapped.1.fastq -fq2 SS6005044_chr22_mapped.2.fastq 
#    X bamToFastq -i SS6005044_chr22_unmapped.sort.bam -fq SS6005044_chr22_unmapped.1.fastq -fq2 SS6005044_chr22_unmapped.2.fastq
samtools fastq -@ 8 SS6005044_chr22_mapped.sort.bam -1 SS6005044_chr22_mapped.1.fastq -2 SS6005044_chr22_mapped.2.fastq -0 /dev/null -s /dev/null -n
samtools fastq -@ 8 SS6005044_chr22_unmapped.sort.bam -1 SS6005044_chr22_unmapped.1.fastq -2 SS6005044_chr22_unmapped.2.fastq -0 /dev/null -s /dev/null -n
# combine both the first and paired reads together from the mapped and unmapped files
cat SS6005044_chr22_mapped.1.fastq SS6005044_chr22_unmapped.1.fastq > SS6005044_chr22.1.fastq
cat SS6005044_chr22_mapped.2.fastq SS6005044_chr22_unmapped.2.fastq > SS6005044_chr22.2.fastq

#BWA
bwa mem /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa SS6005044_chr22_mapped.1.fastq SS6005044_chr22_mapped.2.fastq > aln-pe.sam
bwa aln /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa SS6005044_chr22_mapped.1.fastq > SS6005044_chr22_mapped.1.sai && bwa aln /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa SS6005044_chr22_mapped.2.fastq > SS6005044_chr22_mapped.2.sai && bwa sampe /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa SS6005044_chr22_mapped.1.sai SS6005044_chr22_mapped.2.sai SS6005044_chr22_mapped.1.fastq SS6005044_chr22_mapped.2.fastq > SS6005044_chr22_mapped_1_2.bwa.sam
#BOWTIE2
bowtie2-build /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa chr22
bowtie2 -p 4 -x chr22 -1 SS6005044_chr22_mapped.1.fastq -2 SS6005044_chr22_mapped.2.fastq -S bowtie2_aln_pe.sam





samtools view -Sb aln_pe.sam > aln_pe.bam
samtools sort aln_pe.bam -o aln_pe_sort.bam
samtools index aln_pe_sort.bam
# Adding Read Group tags and indexing bam files
java -jar ${EBROOTPICARD}/picard.jar  AddOrReplaceReadGroups INPUT= aln_pe_sort.bam OUTPUT= aln_pe.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT
samtools index aln_pe.RG.bam
# Marking and removing duplicates
java -jar ${EBROOTPICARD}/picard.jar  MarkDuplicates I= aln_pe.RG.bam O= aln_pe.DR.bam M=aln_pe_output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> aln_pe_logFile.log
samtools index aln_pe.DR.bam


samtools view -Sb SS6005044_chr22_mapped_1_2.bwa.sam > SS6005044_chr22_mapped_1_2.bwa.bam
samtools sort SS6005044_chr22_mapped_1_2.bwa.bam -o SS6005044_chr22_mapped_1_2.bwa.sort.bam
samtools index SS6005044_chr22_mapped_1_2.bwa.sort.bam
# Adding Read Group tags and indexing bam files
java -jar ${EBROOTPICARD}/picard.jar  AddOrReplaceReadGroups INPUT= SS6005044_chr22_mapped_1_2.bwa.sort.bam OUTPUT= SS6005044_chr22_mapped_1_2.bwa.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT
samtools index SS6005044_chr22_mapped_1_2.bwa.RG.bam
# Marking and removing duplicates
java -jar ${EBROOTPICARD}/picard.jar  MarkDuplicates I= SS6005044_chr22_mapped_1_2.bwa.RG.bam O= SS6005044_chr22_mapped_1_2.bwa.DR.bam M=SS6005044_chr22_mapped_1_2.bwa_output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> SS6005044_chr22_mapped_1_2.bwa_logFile.log
samtools index SS6005044_chr22_mapped_1_2.bwa.DR.bam

samtools view -Sb bowtie2_aln_pe.sam > bowtie2_aln_pe.bam
samtools sort bowtie2_aln_pe.bam -o bowtie2_aln_pe_sort.bam
samtools index bowtie2_aln_pe_sort.bam
# Adding Read Group tags and indexing bam files
java -jar ${EBROOTPICARD}/picard.jar  AddOrReplaceReadGroups INPUT= bowtie2_aln_pe_sort.bam OUTPUT= bowtie2_aln_pe.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT
samtools index bowtie2_aln_pe.RG.bam
# Marking and removing duplicates
java -jar ${EBROOTPICARD}/picard.jar  MarkDuplicates I= bowtie2_aln_pe.RG.bam O= bowtie2_aln_pe.DR.bam M=bowtie2_aln_pe_output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> bowtie2_aln_pe_logFile.log
samtools index bowtie2_aln_pe.DR.bam























samtools view -h SS6005044.sorted.bam chr22 > SS6005044_TEST.chr22.sam
samtools view -bS SS6005044_TEST.chr22.sam > SS6005044_TEST.chr22.bam
samtools view -H SS6005044_TEST.chr22.bam
java -jar ${EBROOTPICARD}/picard.jar AddOrReplaceReadGroups I=SS6005044_TEST.chr22.bam O=SS6005044_TEST_add.chr22.sort.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
samtools view -H SS6005044_TEST_add.chr22.sort.bam
        (stage) [umcg-aewijk@gearshift EGAD00001000292]$ samtools view -H SS6005044_TEST_add.chr22.sort.bam
        @HD     VN:1.6  SO:coordinate
        @SQ     SN:chr1 LN:249250621
        @SQ     SN:chr2 LN:243199373
        @SQ     SN:chr3 LN:198022430
        @SQ     SN:chr4 LN:191154276
        @SQ     SN:chr5 LN:180915260
        @SQ     SN:chr6 LN:171115067
        @SQ     SN:chr7 LN:159138663
        @SQ     SN:chrX LN:155270560
        @SQ     SN:chr8 LN:146364022
        @SQ     SN:chr9 LN:141213431
        @SQ     SN:chr10        LN:135534747
        @SQ     SN:chr11        LN:135006516
        @SQ     SN:chr12        LN:133851895
        @SQ     SN:chr13        LN:115169878
        @SQ     SN:chr14        LN:107349540
        @SQ     SN:chr15        LN:102531392
        @SQ     SN:chr16        LN:90354753
        @SQ     SN:chr17        LN:81195210
        @SQ     SN:chr18        LN:78077248
        @SQ     SN:chr20        LN:63025520
        @SQ     SN:chr19        LN:59128983
        @SQ     SN:chr22        LN:51304566
        @SQ     SN:chr21        LN:48129895
        @SQ     SN:chrM LN:16571
        @RG     ID:4    LB:lib1 PL:illumina     SM:20   PU:unit1
        @PG     ID:CASAVA       VN:CASAVA-1.9.0a1_110909        CL:/illumina/development/casava/VariantCalling_Version2_5a/bin/configureBuild.pl --targets all bam --inSampleDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Aligned/D0LCYACXX/Sample_SS6005044 --inSampleDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Aligned/D0LR3ACXX/Sample_SS6005044 --outDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Assembly --samtoolsRefFile=/isilon/Genomes/FASTA_UCSC/HumanNCBI37_UCSC/HumanNCBI37_UCSC_XX.fa --indelsSaveTempFiles --jobsLimit=80 --variantsPrintUsedAlleleCounts --variantsWriteRealigned --sortKeepAllReads --bamChangeChromLabels=OFF --sgeQueue=all.q

#http://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation
java -jar ${EBROOTPICARD}/picard.jar FixMateInformation I=SS6005044_TEST_add.chr22.sort.bam O=SS6005044_TEST_FIXM.chr22.sort.bam 
samtools view -H SS6005044_TEST_FIXM.chr22.sort.bam

samtools view -H SS6005044_TEST_FIXM.chr22.sort.bam | grep @HD > new_header.txt
samtools view -H SS6005044_TEST_FIXM.chr22.sort.bam | grep @PG >> new_header.txt
samtools view -H SS6005044_TEST_FIXM.chr22.sort.bam | grep @RG >> new_header.txt
samtools view -H SS6005044_TEST_FIXM.chr22.sort.bam | grep chr22  >> new_header.txt

samtools reheader new_header.txt SS6005044_TEST_FIXM.chr22.sort.bam > SS6005044_TEST_FIXM.chr22.reheadered.bam
samtools index SS6005044_TEST_FIXM.chr22.reheadered.bam
samtools view -H SS6005044_TEST_FIXM.chr22.reheadered.bam
        @HD     VN:1.6  SO:coordinate
        @PG     ID:CASAVA       VN:CASAVA-1.9.0a1_110909        CL:/illumina/development/casava/VariantCalling_Version2_5a/bin/configureBuild.pl --targets all bam --inSampleDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Aligned/D0LCYACXX/Sample_SS6005044 --inSampleDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Aligned/D0LR3ACXX/Sample_SS6005044 --outDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Assembly --samtoolsRefFile=/isilon/Genomes/FASTA_UCSC/HumanNCBI37_UCSC/HumanNCBI37_UCSC_XX.fa --indelsSaveTempFiles --jobsLimit=80 --variantsPrintUsedAlleleCounts --variantsWriteRealigned --sortKeepAllReads --bamChangeChromLabels=OFF --sgeQueue=all.q
        @RG     ID:4    LB:lib1 PL:illumina     SM:20   PU:unit1
        @SQ     SN:chr22        LN:51304566
        @PG     ID:samtools     PN:samtools     PP:CASAVA       VN:1.9  CL:samtools reheader new_header.txt SS6005044_TEST_FIXM.chr22.sort.bam


gatk Mutect2 -R /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa -I SS6005044_TEST_FIXM.chr22.reheadered.bam -O SS6005044_TEST_FIXM.chr22.reheadered_somatic.vcf.gz
gatk FilterMutectCalls -R /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa -V SS6005044_TEST_FIXM.chr22.reheadered_somatic.vcf.gz -O SS6005044_TEST_FIXM.chr22.reheadered_somatic.vcf.gz





















#!/usr/bin/bash

cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/

ml Anaconda3/5.3.0
ml SAMtools/1.9-foss-2018b
ml GATK/4.1.4.1-Java-8-LTS
ml BCFtools/1.11-GCCcore-7.3.0
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS

CHR=chr21
FILE=SS6005044


#Preparation of the reference sequence
bwa index -a is /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHR}.fa
# Create a reference dictionary
java -jar ${EBROOTPICARD}/picard.jar CreateSequenceDictionary R= /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHR}.fa O= /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHR}.dict
# Index the reference sequence with Samtools
samtools faidx /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHR}.fa



samtools view -h ${FILE}.sorted.bam ${CHR} > ${FILE}_TEST.${CHR}.sam
samtools view -bS ${FILE}_TEST.${CHR}.sam > ${FILE}_TEST.${CHR}.bam
samtools view -H ${FILE}_TEST.${CHR}.bam
java -jar ${EBROOTPICARD}/picard.jar AddOrReplaceReadGroups I=${FILE}_TEST.${CHR}.bam O=${FILE}_TEST_add.${CHR}.sort.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
samtools view -H ${FILE}_TEST_add.${CHR}.sort.bam
#        (stage) [umcg-aewijk@gearshift EGAD00001000292]$ samtools view -H SS6005044_TEST_add.chr22.sort.bam
#        @HD     VN:1.6  SO:coordinate
#        @SQ     SN:chr1 LN:249250621
#        @SQ     SN:chr2 LN:243199373
#        @SQ     SN:chr3 LN:198022430
#        @SQ     SN:chr4 LN:191154276
#        @SQ     SN:chr5 LN:180915260
#        @SQ     SN:chr6 LN:171115067
#        @SQ     SN:chr7 LN:159138663
#        @SQ     SN:chrX LN:155270560
#        @SQ     SN:chr8 LN:146364022
#        @SQ     SN:chr9 LN:141213431
#        @SQ     SN:chr10        LN:135534747
#        @SQ     SN:chr11        LN:135006516
#        @SQ     SN:chr12        LN:133851895
#        @SQ     SN:chr13        LN:115169878
#        @SQ     SN:chr14        LN:107349540
#        @SQ     SN:chr15        LN:102531392
#        @SQ     SN:chr16        LN:90354753
#        @SQ     SN:chr17        LN:81195210
#        @SQ     SN:chr18        LN:78077248
#        @SQ     SN:chr20        LN:63025520
#        @SQ     SN:chr19        LN:59128983
#        @SQ     SN:chr22        LN:51304566
#        @SQ     SN:chr21        LN:48129895
#        @SQ     SN:chrM LN:16571
#        @RG     ID:4    LB:lib1 PL:illumina     SM:20   PU:unit1
#        @PG     ID:CASAVA       VN:CASAVA-1.9.0a1_110909        CL:/illumina/development/casava/VariantCalling_Version2_5a/bin/configureBuild.pl --targets all bam --inSampleDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Aligned/D0LCYACXX/Sample_SS6005044 --inSampleDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Aligned/D0LR3ACXX/Sample_SS6005044 --outDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Assembly --samtoolsRefFile=/isilon/Genomes/FASTA_UCSC/HumanNCBI37_UCSC/HumanNCBI37_UCSC_XX.fa --indelsSaveTempFiles --jobsLimit=80 --variantsPrintUsedAlleleCounts --variantsWriteRealigned --sortKeepAllReads --bamChangeChromLabels=OFF --sgeQueue=all.q

#http://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation
java -jar ${EBROOTPICARD}/picard.jar FixMateInformation I=${FILE}_TEST_add.${CHR}.sort.bam O=${FILE}_TEST_FIXM.${CHR}.sort.bam 
samtools view -H ${FILE}_TEST_FIXM.${CHR}.sort.bam

samtools view -H ${FILE}_TEST_FIXM.${CHR}.sort.bam | grep @HD > new_header.txt
samtools view -H ${FILE}_TEST_FIXM.${CHR}.sort.bam | grep @PG >> new_header.txt
samtools view -H ${FILE}_TEST_FIXM.${CHR}.sort.bam | grep @RG >> new_header.txt
samtools view -H ${FILE}_TEST_FIXM.${CHR}.sort.bam | grep ${CHR}  >> new_header.txt

samtools reheader new_header.txt ${FILE}_TEST_FIXM.${CHR}.sort.bam > ${FILE}_TEST_FIXM.${CHR}.reheadered.bam
samtools index ${FILE}_TEST_FIXM.${CHR}.reheadered.bam
samtools view -H ${FILE}_TEST_FIXM.${CHR}.reheadered.bam
#        @HD     VN:1.6  SO:coordinate
#        @PG     ID:CASAVA       VN:CASAVA-1.9.0a1_110909        CL:/illumina/development/casava/VariantCalling_Version2_5a/bin/configureBuild.pl --targets all bam --inSampleDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Aligned/D0LCYACXX/Sample_SS6005044 --inSampleDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Aligned/D0LR3ACXX/Sample_SS6005044 --outDir=/isilon/RUO/Projects/CRUK_Fitzgibbon/SS6005044/Assembly --samtoolsRefFile=/isilon/Genomes/FASTA_UCSC/HumanNCBI37_UCSC/HumanNCBI37_UCSC_XX.fa --indelsSaveTempFiles --jobsLimit=80 --variantsPrintUsedAlleleCounts --variantsWriteRealigned --sortKeepAllReads --bamChangeChromLabels=OFF --sgeQueue=all.q
#        @RG     ID:4    LB:lib1 PL:illumina     SM:20   PU:unit1
#        @SQ     SN:chr22        LN:51304566
#        @PG     ID:samtools     PN:samtools     PP:CASAVA       VN:1.9  CL:samtools reheader new_header.txt SS6005044_TEST_FIXM.chr22.sort.bam


gatk Mutect2 -R /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHR}.fa -I ${FILE}_TEST_FIXM.${CHR}.reheadered.bam -O ${FILE}_TEST_FIXM.${CHR}.reheadered_somatic.vcf.gz
gatk FilterMutectCalls -R /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHR}.fa -V ${FILE}_TEST_FIXM.${CHR}.reheadered_somatic.vcf.gz -O ${FILE}_TEST_FIXM.${CHR}.reheadered_somatic.vcf.gz








https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/01_alignment.html




