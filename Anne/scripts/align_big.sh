#!/usr/bin/bash

#SBATCH --job-name=align_big
#SBATCH --output=align_big.out
#SBATCH --error=align_big.err
#SBATCH --time=167:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

NUMBER=5044
FILE_NUM=SS600${NUMBER}
GENOOM=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/chrall.fa
PATH=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/

ml SAMtools/1.9-foss-2018b
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS

echo "BEGIN"
echo "${NUMBER}"

mkdir /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${NUMBER}

#cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292
#samtools view -h ${FILE_NUM}.sorted.bam ${CHROM} > ${CHROM}/${NUMBER}/${FILE_NUM}.${CHROM}.sam



mkdir /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${NUMBER}/bowtie
mkdir /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${NUMBER}/bwa_aln
mkdir /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${NUMBER}/bwa_mem
mkdir /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${NUMBER}/QC

#samtools view -bS ${FILE_NUM}.${CHROM}.sam > ${FILE_NUM}.${CHROM}.bam
#samtools sort ${FILE_NUM}.${CHROM}.bam -o ${FILE_NUM}_${CHROM}.sort.bam
#samtools sort -n ${FILE_NUM}.${CHROM}.bam -o ${FILE_NUM}_name_${CHROM}.sort.bam


#SS6005044.sorted.bam
#reads that mapped properly as pairs
samtools view -u -f 1 -F 12 ${PATH}${FILE_NUM}.sorted.bam > ${PATH}${NUMBER}/${FILE_NUM}_map_map.bam
# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 ${PATH}${FILE_NUM}.sorted.bam > ${PATH}${NUMBER}/${FILE_NUM}_unmap_map.bam
# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 ${PATH}${FILE_NUM}.sorted.bam > ${PATH}${NUMBER}/${FILE_NUM}_map_unmap.bam
# R1 & R2 unmapped
samtools view -u -f 12 -F 256 ${PATH}${FILE_NUM}.sorted.bam > ${PATH}${NUMBER}/${FILE_NUM}_unmap_unmap.bam

cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/${NUMBER}
# merge the three files
samtools merge -u ${PATH}${NUMBER}${FILE_NUM}_unmapped.bam ${PATH}${NUMBER}${FILE_NUM}_unmap_map.bam ${PATH}${NUMBER}${FILE_NUM}_map_unmap.bam ${PATH}${NUMBER}${FILE_NUM}_unmap_unmap.bam
# BAM files must be resorted so that they are ordered by read ID instead of location in the reference
samtools sort -n ${PATH}${NUMBER}${FILE_NUM}_map_map.bam -o ${PATH}${NUMBER}${FILE_NUM}_mapped.sort.bam
samtools sort -n ${PATH}${NUMBER}${FILE_NUM}_unmapped.bam -o ${PATH}${NUMBER}${FILE_NUM}_unmapped.sort.bam
# it is a good idea to check that you have the correct number of reads and no redundancy. You can summarize the original BAM file to get an idea of where you started.
samtools flagstat ${PATH}${NUMBER}${FILE_NUM}.sorted.bam
samtools view -c ${PATH}${NUMBER}${FILE_NUM}_mapped.sort.bam
samtools view -c ${PATH}${NUMBER}${FILE_NUM}_unmapped.sort.bam
# extract the FASTQ reads into two paired read files
#    X bamToFastq -i SS6005044_chr22_mapped.sort.bam -fq SS6005044_chr22_mapped.1.fastq -fq2 SS6005044_chr22_mapped.2.fastq 
#    X bamToFastq -i SS6005044_chr22_unmapped.sort.bam -fq SS6005044_chr22_unmapped.1.fastq -fq2 SS6005044_chr22_unmapped.2.fastq
samtools fastq -@ 8 ${PATH}${NUMBER}${FILE_NUM}_mapped.sort.bam -1 ${PATH}${NUMBER}${FILE_NUM}_mapped.1.fastq -2 ${PATH}${NUMBER}${FILE_NUM}_mapped.2.fastq -0 /dev/null -s /dev/null -n
samtools fastq -@ 8 ${PATH}${NUMBER}${FILE_NUM}_unmapped.sort.bam -1 ${PATH}${NUMBER}${FILE_NUM}_unmapped.1.fastq -2 ${PATH}${NUMBER}${FILE_NUM}_unmapped.2.fastq -0 /dev/null -s /dev/null -n
# combine both the first and paired reads together from the mapped and unmapped files
cat ${PATH}${NUMBER}${FILE_NUM}_mapped.1.fastq ${PATH}${NUMBER}${FILE_NUM}_unmapped.1.fastq > ${PATH}${NUMBER}${FILE_NUM}.1.fastq
cat ${PATH}${NUMBER}${FILE_NUM}_mapped.2.fastq ${PATH}${NUMBER}${FILE_NUM}_unmapped.2.fastq > ${PATH}${NUMBER}${FILE_NUM}.2.fastq

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
bwa mem ${GENOOM} ${PATH}${NUMBER}${FILE_NUM}_mapped.1.fastq ${PATH}${NUMBER}${FILE_NUM}_mapped.2.fastq > ${PATH}${NUMBER}/bwa_mem/mem_${FILE_NUM}.sam
align_test ${PATH}${NUMBER}/bwa_mem/mem_${FILE_NUM}

bwa aln ${GENOOM} ${PATH}${NUMBER}${FILE_NUM}_mapped.1.fastq > ${PATH}${NUMBER}${FILE_NUM}_mapped.1.sai && bwa aln ${GENOOM} ${PATH}${NUMBER}${FILE_NUM}_mapped.2.fastq > ${PATH}${NUMBER}${FILE_NUM}_mapped.2.sai && bwa sampe ${GENOOM} ${PATH}${NUMBER}${FILE_NUM}_mapped.1.sai ${PATH}${NUMBER}${FILE_NUM}_mapped.2.sai ${PATH}${NUMBER}${FILE_NUM}_mapped.1.fastq ${PATH}${NUMBER}${FILE_NUM}_mapped.2.fastq > ${PATH}${NUMBER}/bwa_aln/aln_${FILE_NUM}.sam
align_test ${PATH}${NUMBER}/bwa_aln/aln_${FILE_NUM}

#BOWTIE2
bowtie2-build ${GENOOM} ${PATH}${NUMBER}/bowtie/GENOOM
bowtie2 -p 4 -x ${PATH}${NUMBER}/bowtie/GENOOM -1 ${PATH}${NUMBER}${FILE_NUM}_mapped.1.fastq -2 ${PATH}${NUMBER}${FILE_NUM}_mapped.2.fastq -S ${PATH}${NUMBER}/bowtie/bowtie2_${FILE_NUM}.sam
align_test ${PATH}${NUMBER}/bowtie/bowtie2_${FILE_NUM}


echo "EIND"