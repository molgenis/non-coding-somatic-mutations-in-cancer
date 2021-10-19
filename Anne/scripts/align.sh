#!/usr/bin/bash

NUMBER=5044
CHROM=chr22

ml SAMtools/1.9-foss-2018b
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS

echo "BEGIN"
echo "${NUMBER}"

mkdir /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${CHROM}/${NUMBER}

cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292
samtools view -h SS600${NUMBER}.sorted.bam ${CHROM} > ${CHROM}/${NUMBER}/SS600${NUMBER}.${CHROM}.sam


cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${CHROM}/${NUMBER}
mkdir /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${CHROM}/${NUMBER}/bowtie
mkdir /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${CHROM}/${NUMBER}/bwa_aln
mkdir /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${CHROM}/${NUMBER}/bwa_mem
mkdir /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${CHROM}/${NUMBER}/QC

samtools view -bS SS600${NUMBER}.${CHROM}.sam > SS600${NUMBER}.${CHROM}.bam
samtools sort SS600${NUMBER}.${CHROM}.bam -o SS600${NUMBER}_${CHROM}.sort.bam
samtools sort -n SS600${NUMBER}.${CHROM}.bam -o SS600${NUMBER}_name_${CHROM}.sort.bam


#reads that mapped properly as pairs
samtools view -u -f 1 -F 12 SS600${NUMBER}_${CHROM}.sort.bam > SS600${NUMBER}_${CHROM}_map_map.bam
# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 SS600${NUMBER}_${CHROM}.sort.bam > SS600${NUMBER}_${CHROM}_unmap_map.bam
# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 SS600${NUMBER}_${CHROM}.sort.bam > SS600${NUMBER}_${CHROM}_map_unmap.bam
# R1 & R2 unmapped
samtools view -u -f 12 -F 256 SS600${NUMBER}_${CHROM}.sort.bam > SS600${NUMBER}_${CHROM}_unmap_unmap.bam
# merge the three files
samtools merge -u SS600${NUMBER}_${CHROM}_unmapped.bam SS600${NUMBER}_${CHROM}_unmap_map.bam SS600${NUMBER}_${CHROM}_map_unmap.bam SS600${NUMBER}_${CHROM}_unmap_unmap.bam
# BAM files must be resorted so that they are ordered by read ID instead of location in the reference
samtools sort -n SS600${NUMBER}_${CHROM}_map_map.bam -o SS600${NUMBER}_${CHROM}_mapped.sort.bam
samtools sort -n SS600${NUMBER}_${CHROM}_unmapped.bam -o SS600${NUMBER}_${CHROM}_unmapped.sort.bam
# it is a good idea to check that you have the correct number of reads and no redundancy. You can summarize the original BAM file to get an idea of where you started.
samtools flagstat SS600${NUMBER}_${CHROM}.sort.bam
samtools view -c SS600${NUMBER}_${CHROM}_mapped.sort.bam
samtools view -c SS600${NUMBER}_${CHROM}_unmapped.sort.bam
# extract the FASTQ reads into two paired read files
#    X bamToFastq -i SS6005044_chr22_mapped.sort.bam -fq SS6005044_chr22_mapped.1.fastq -fq2 SS6005044_chr22_mapped.2.fastq 
#    X bamToFastq -i SS6005044_chr22_unmapped.sort.bam -fq SS6005044_chr22_unmapped.1.fastq -fq2 SS6005044_chr22_unmapped.2.fastq
samtools fastq -@ 8 SS600${NUMBER}_${CHROM}_mapped.sort.bam -1 SS600${NUMBER}_${CHROM}_mapped.1.fastq -2 SS600${NUMBER}_${CHROM}_mapped.2.fastq -0 /dev/null -s /dev/null -n
samtools fastq -@ 8 SS600${NUMBER}_${CHROM}_unmapped.sort.bam -1 SS600${NUMBER}_${CHROM}_unmapped.1.fastq -2 SS600${NUMBER}_${CHROM}_unmapped.2.fastq -0 /dev/null -s /dev/null -n
# combine both the first and paired reads together from the mapped and unmapped files
cat SS600${NUMBER}_${CHROM}_mapped.1.fastq SS600${NUMBER}_${CHROM}_unmapped.1.fastq > SS600${NUMBER}_${CHROM}.1.fastq
cat SS600${NUMBER}_${CHROM}_mapped.2.fastq SS600${NUMBER}_${CHROM}_unmapped.2.fastq > SS600${NUMBER}_${CHROM}.2.fastq

align_test() {
    samtools view -Sb ${1}.sam > ${1}.bam
    samtools sort ${1}.bam -o ${1}_sort.bam
    samtools index ${1}_sort.bam
    # Adding Read Group tags and indexing bam files
    java -jar ${EBROOTPICARD}/picard.jar  AddOrReplaceReadGroups INPUT= ${1}_sort.bam OUTPUT= ${1}.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT
    samtools index ${1}.RG.bam
    # Marking and removing duplicates
    java -jar ${EBROOTPICARD}/picard.jar  MarkDuplicates I= ${1}.RG.bam O= ${1}.DR.bam M=${1}_output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> ${1}_logFile.log
    samtools index ${1}.DR.bam
}

#BWA = Burrows-Wheeler Aligner
#https://ucdavis-bioinformatics-training.github.io/2017-August-Variant-Analysis-Workshop/wednesday/alignment.html
bwa mem /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHROM}.fa SS600${NUMBER}_${CHROM}_mapped.1.fastq SS600${NUMBER}_${CHROM}_mapped.2.fastq > bwa_mem/aln_pe.sam
align_test bwa_mem/aln_pe

bwa aln /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHROM}.fa SS600${NUMBER}_${CHROM}_mapped.1.fastq > SS600${NUMBER}_${CHROM}_mapped.1.sai && bwa aln /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHROM}.fa SS600${NUMBER}_${CHROM}_mapped.2.fastq > SS600${NUMBER}_${CHROM}_mapped.2.sai && bwa sampe /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHROM}.fa SS600${NUMBER}_${CHROM}_mapped.1.sai SS600${NUMBER}_${CHROM}_mapped.2.sai SS600${NUMBER}_${CHROM}_mapped.1.fastq SS600${NUMBER}_${CHROM}_mapped.2.fastq > bwa_aln/bwa_aln_${NUMBER}_${CHROM}_mapped_1_2.sam
align_test bwa_aln/bwa_aln_${NUMBER}_${CHROM}_mapped_1_2

#BOWTIE2
bowtie2-build /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/${CHROM}.fa ${CHROM}
bowtie2 -p 4 -x ${CHROM} -1 SS600${NUMBER}_${CHROM}_mapped.1.fastq -2 SS600${NUMBER}_${CHROM}_mapped.2.fastq -S bowtie/bowtie2_aln_pe.sam
align_test bowtie/bowtie2_aln_pe


echo "EIND"



