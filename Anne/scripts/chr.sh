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
FILE=SS6005043


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