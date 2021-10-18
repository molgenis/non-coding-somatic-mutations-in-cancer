#!/usr/bin/bash

#FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/SS6004099.sorted.bam
#TYPE=tumor #tumor or normal

TUMOR=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/
NORMAL=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/

PATH_REF=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/
REF=chr22.fa

CHR_DICT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.dict
CHR=chr22

VCF_PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/
VCF_NAME=somatic-b37_Mutect2-WGS-panel-b37.vcf
PON=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/merge_somatic-b37_Mutect2-WGS-panel-b37.vcf

ml Anaconda3/5.3.0
ml SAMtools/1.9-foss-2018b
ml GATK/4.1.4.1-Java-8-LTS
ml BCFtools/1.11-GCCcore-7.3.0
ml BWA/0.7.17-GCCcore-7.3.0
ml Bowtie2/2.3.4.2-foss-2018b
ml picard/2.20.5-Java-11-LTS


sort_index() {
    samtools sort ${FILE}.bam -o ${FILE}.sorted.bam
    samtools index ${FILE}.sorted.bam
}

filter_on_chr(){
    samtools index ${1}.sorted.bam
    samtools view -h ${1}.sorted.bam ${CHR} > ${1}.${CHR}.sam
    samtools view -bS ${1}.${CHR}.sam > ${1}.${CHR}.bam
    printf "filter check: %s" "${1}"
    samtools sort ${1}.${CHR}.bam -o ${1}.${CHR}.sorted.bam
    samtools index ${1}.${CHR}.sorted.bam
}

align_again() {
    #bam omzetten naar  fastq
    samtools fastq ${1}.${CHR}.bam > ${1}.${CHR}.fastq
    #fastq zippen
    gzip ${1}.${CHR}.fastq
    #Preparation of the reference sequence
    bwa index -a is ${PATH_REF}${REF}
    #Create a reference dictionary
    java -jar ${EBROOTPICARD}/picard.jar CreateSequenceDictionary R= ${PATH_REF}${REF}  O= ${CHR_DICT}
    #Index the reference sequence with Samtools
    samtools faidx ${PATH_REF}${REF}
    #Alignment of pre-processed reads to the reference genome with BWA aln
    bwa aln ${PATH_REF}${REF} ${1}.${CHR}.fastq.gz -n 0.1 -l 1000 > ${1}.${CHR}.sai
    bwa samse ${PATH_REF}${REF} ${1}.${CHR}.sai ${1}.${CHR}.fastq.gz -f ${1}.${CHR}.sam
    #Converting sam file to bam file
    samtools view -Sb ${1}.${CHR}.sam > ${1}.${CHR}.bam
    samtools view ${1}.${CHR}.bam | less -S
    samtools view -H ${1}.${CHR}.bam
    #Sorting and indexing the bam file
    samtools sort ${1}.${CHR}.bam -o ${1}.${CHR}.sort.bam
    samtools index ${1}.${CHR}.sort.bam
    #Adding Read Group tags and indexing bam files
    java -jar ${EBROOTPICARD}/picard.jar  AddOrReplaceReadGroups INPUT= ${1}.${CHR}.sort.bam OUTPUT= ${1}.${CHR}.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT
    samtools index ${1}.${CHR}.RG.bam
    #Marking and removing duplicates
    java -jar ${EBROOTPICARD}/picard.jar  MarkDuplicates I= ${1}.${CHR}.RG.bam O= ${1}.${CHR}.DR.bam M=output_metrics_${1}_${CHR}.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> logFile_${1}_${CHR}.log
    samtools index ${1}.${CHR}.DR.bam
}

change_reference_genome() {
    # http://seqanswers.com/forums/showthread.php?t=22504
    sed -e 's/##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/' -e 's/>MT/>chrM/' ${PATH_REF}${REF} > ${PATH_REF}New${REF}
    #awk '/^>/ {P=index($0,"chrM")==0} {if(P) print} ' NEW_genome.fa > without_chrM_genome.fa
}

change_VCF() {
    # header
    bcftools view --header-only ${VCF_PON}${VCF_NAME} | sed 's/##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/' | sed 's/##contig=<ID=MT/##contig=<ID=chrM/' > ${VCF_PON}head_${VCF_NAME}
    # no headers
    bcftools view --header-only ${VCF_PON}${VCF_NAME} | sed 's/##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/' | sed 's/##contig=<ID=MT/##contig=<ID=chrM/' > ${VCF_PON}nohead_${VCF_NAME}
    # header en no header plakken
    cat ${VCF_PON}head_${VCF_NAME} ${VCF_PON}nohead_${VCF_NAME} >> ${VCF_PON}merge_${VCF_NAME}
    # index
    gatk IndexFeatureFile -I ${VCF_PON}merge_${VCF_NAME}
}

change_sample_names() {
    #for f in *.DR.bam; do echo -ne "$f\t" ; samtools view -H $f | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq; done
    for f in *.DR.bam
    do 
        echo -ne "$f\t" 
        samtools view -H $f  | sed "s/SM:[^\t]*/SM:${f}/g" | samtools reheader - $f > SN_$f
    done
}

change_sample_name() {
    samtools view -H ${1}.${CHR}.DR.bam  | sed "s/SM:[^\t]*/SM:${1}/g" | samtools reheader - ${1}.${CHR}.DR.bam > SN_${1}.${CHR}.DR.bam
    samtools index SN_${1}.${CHR}.DR.bam
}

Mutect2_one_sample(){
    # standard
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${CHR}.DR.bam -O unfiltered_${1}_${CHR}_somatic.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${CHR}_somatic.vcf.gz -O filtered_${1}_${CHR}_somatic.vcf.gz
    # sample name
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${CHR}.DR.bam -${2} ${1}.${CHR}.DR.bam -O unfiltered_${1}_${CHR}_somatic_name.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${CHR}_somatic_name.vcf.gz -O filtered_${1}_${CHR}_somatic_name.vcf.gz
    # sample name and PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${CHR}.DR.bam -${2} ${1}.${CHR}.DR.bam -O unfiltered_${1}_${CHR}_somatic_name_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${CHR}_somatic_name_PON.vcf.gz -O filtered_${1}_${CHR}_somatic_name_PON.vcf.gz
    # PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${CHR}.DR.bam -O unfiltered_${1}_${CHR}_somatic_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${CHR}_somatic_PON.vcf.gz -O filtered_${1}_${CHR}_somatic_PON.vcf.gz   
}

Mutect2_two_samples() {
    # standard
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${CHR}.DR.bam -I SN_${2}.${CHR}.DR.bam -O unfiltered_${1}_${2}_${CHR}_somatic.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${2}_${CHR}_somatic.vcf.gz -O filtered_${1}_${2}_${CHR}_somatic.vcf.gz
    # sample name
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${CHR}.DR.bam -tumor ${1}.${CHR}.DR.bam -I SN_${2}.${CHR}.DR.bam -normal ${2}.${CHR}.DR.bam -O unfiltered_${1}_${2}_${CHR}_somatic_name.vcf.gz 
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${2}_${CHR}_somatic_name.vcf.gz -O filtered_${1}_${2}_${CHR}_somatic_name.vcf.gz
    # sample name and PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${CHR}.DR.bam -tumor ${1}.${CHR}.DR.bam -I SN_${2}.${CHR}.DR.bam -normal ${2}.${CHR}.DR.bam -O unfiltered_${1}_${2}_${CHR}_somatic_name_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${2}_${CHR}_somatic_name_PON.vcf.gz -O filtered_${1}_${2}_${CHR}_somatic_name_PON.vcf.gz
    # PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_${1}.${CHR}.DR.bam -I SN_${2}.${CHR}.DR.bam -O unfiltered_${1}_${2}_${CHR}_somatic_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_${1}_${2}_${CHR}_somatic_PON.vcf.gz -O filtered_${1}_${2}_${CHR}_somatic_PON.vcf.gz
}

Mutect2_three_samples(){
    # standard
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_$1.${CHR}.DR.bam -I SN_$2.${CHR}.DR.bam -I SN_$3.${CHR}.DR.bam -O unfiltered_$1_$2_$3_${CHR}_somatic.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_$1_$2_$3_${CHR}_somatic.vcf.gz -O filtered_$1_$2_$3_${CHR}_somatic.vcf.gz
    # sample name
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_$1.${CHR}.DR.bam -I SN_$2.${CHR}.DR.bam -I SN_$3.${CHR}.DR.bam -normal $3.${CHR}.DR.bam -O unfiltered_$1_$2_$3_${CHR}_somatic_name.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_$1_$2_$3_${CHR}_somatic_name.vcf.gz -O filtered_$1_$2_$3_${CHR}_somatic_name.vcf.gz
    # sample name and PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_$1.${CHR}.DR.bam -I SN_$2.${CHR}.DR.bam -I SN_$3.${CHR}.DR.bam -normal $3.${CHR}.DR.bam -O unfiltered_$1_$2_$3_${CHR}_somatic_name_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_$1_$2_$3_${CHR}_somatic_name_PON.vcf.gz -O filtered_$1_$2_$3_${CHR}_somatic_name_PON.vcf.gz
    # PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SN_$1.${CHR}.DR.bam -I SN_$2.${CHR}.DR.bam -I SN_$3.${CHR}.DR.bam -O unfiltered_$1_$2_$3_${CHR}_somatic_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V unfiltered_$1_$2_$3_${CHR}_somatic_PON.vcf.gz -O filtered_$1_$2_$3_${CHR}_somatic_PON.vcf.gz
}

#Filter_with_hand() {
#}

cd /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/

#printf "####change_sample_name\n"
#printf "###########5044\n"
#change_sample_name SS6005044 #1=file
#printf "###########5043\n"
#change_sample_name SS6005043 #1=file
#printf "###########5042\n"
#change_sample_name SS6005042 #1=file

printf "\n####Mutect2_one_sample\n"
#printf "###########5044\n"
#Mutect2_one_sample SS6005044 tumor #1=file 2=type
printf "###########5043\n"
Mutect2_one_sample SS6005043 tumor #1=file 2=type
printf "###########5042\n"
Mutect2_one_sample SS6005042 normal #1=file 2=type

#printf "\n####Mutect2_two_sample\n"
#printf "###########5044-5042\n"
#Mutect2_two_samples SS6005044 SS6005042 #1=tumor 2=normal
#printf "###########5043-5042\n"
#Mutect2_two_samples SS6005043 SS6005042

#printf "\n####Mutect2_three_sample\n"
#printf "###########all\n"
#Mutect2_three_samples SS6005044 SS6005043 SS6005042 #1=tumor 2=tumor 3=normal




