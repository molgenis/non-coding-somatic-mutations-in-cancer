#!/usr/bin/bash

FILE=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/
TYPE=tumor #tumor or normal

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
    samtools view -h ${FILE}.sorted.bam ${CHR} > ${FILE}.${CHR}.sam
    samtools view -bS ${FILE}.${CHR}.sam > ${FILE}.${CHR}.bam
}

align_again() {
    #bam omzetten naar  fastq
    samtools fastq ${FILE}.${CHR}.bam > ${FILE}.${CHR}.fastq
    #fastq zippen
    gzip ${FILE}.${CHR}.fastq
    #Preparation of the reference sequence
    bwa index -a is ${PATH_REF}${REF}
    #Create a reference dictionary
    java -jar ${EBROOTPICARD}/picard.jar CreateSequenceDictionary R= ${PATH_REF}${REF}  O= ${CHR_DICT}
    #Index the reference sequence with Samtools
    samtools faidx ${PATH_REF}${REF}
    #Alignment of pre-processed reads to the reference genome with BWA aln
    bwa aln ${PATH_REF}${REF} ${FILE}.${CHR}.fastq.gz -n 0.1 -l 1000 > ${FILE}.${CHR}.sai
    bwa samse ${PATH_REF}${REF} ${FILE}.${CHR}.sai ${FILE}.${CHR}.fastq.gz -f ${FILE}.${CHR}.sam
    #Converting sam file to bam file
    samtools view -Sb ${FILE}.${CHR}.sam > ${FILE}.${CHR}.bam
    samtools view ${FILE}.${CHR}.bam | less -S
    samtools view -H ${FILE}.${CHR}.bam
    #Sorting and indexing the bam file
    samtools sort ${FILE}.${CHR}.bam -o ${FILE}.${CHR}.sort.bam
    samtools index ${FILE}.${CHR}.sort.bam
    #Adding Read Group tags and indexing bam files
    java -jar ${EBROOTPICARD}/picard.jar  AddOrReplaceReadGroups INPUT= ${FILE}.${CHR}.sort.bam OUTPUT= ${FILE}.${CHR}.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT
    samtools index ${FILE}.${CHR}.RG.bam
    #Marking and removing duplicates
    java -jar ${EBROOTPICARD}/picard.jar  MarkDuplicates I= ${FILE}.${CHR}.RG.bam O= ${FILE}.${CHR}.DR.bam M=output_metrics_${FILE}_${CHR}.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> logFile_${FILE}_${CHR}.log
    samtools index ${FILE}.${CHR}.DR.bam
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

change_sample_name() {
    #for f in *.DR.bam; do echo -ne "$f\t" ; samtools view -H $f | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq; done
    for f in *.DR.bam
    do 
    echo -ne "$f\t" 
    samtools view -H $f  | sed "s/SM:[^\t]*/SM:${f}/g" | samtools reheader - $f > SM_$f
    done
}

Mutect2_one_sample(){
    # standard
    gatk Mutect2 -R ${PATH_REF}${REF} -I SM_${FILE}.${CHR}.DR.bam -O ${FILE}_${CHR}_somatic.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${FILE}_${CHR}_somatic.vcf.gz -O filtered_${FILE}_${CHR}_somatic.vcf.gz
    # sample name
    gatk Mutect2 -R ${PATH_REF}${REF} -I SM_${FILE}.${CHR}.DR.sort.bam -${TYPE} ${FILE}.${CHR}.DR.bam -O ${FILE}_${CHR}_somatic_name.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${FILE}_${CHR}_somatic_name.vcf.gz -O filtered_${FILE}_${CHR}_somatic_name.vcf.gz
    # sample name and PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SM_${FILE}.${CHR}.DR.sort.bam -${TYPE} ${FILE}.${CHR}.DR.bam -O ${FILE}_${CHR}_somatic_name_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${FILE}_${CHR}_somatic_name_PON.vcf.gz -O filtered_${FILE}_${CHR}_somatic_name_PON.vcf.gz
    # PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SM_${FILE}.${CHR}.DR.sort.bam -O ${FILE}_${CHR}_somatic_PON.vcf.gz --panel-of-normals ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${FILE}_${CHR}_somatic_PON.vcf.gz -O filtered_${FILE}_${CHR}_somatic_PON.vcf.gz   
}

Mutect2_two_samples() {
    # standard
    gatk Mutect2 -R ${PATH_REF}${REF} -I SM_${TUMOR}.${CHR}.DR.sort.bam -I SM_${FILE}.${CHR}.DR.sort.bam -O ${TUMOR}_${NORMAL}_${CHR}_somatic.vcf.gz
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${TUMOR}_${NORMAL}_${CHR}_somatic.vcf.gz -O filtered_${TUMOR}_${NORMAL}_${CHR}_somatic.vcf.gz
    # sample name
    gatk Mutect2 -R ${PATH_REF}${REF} -I SM_${TUMOR}.${CHR}.DR.sort.bam -tumor ${TUMOR}.${CHR}.DR.bam -I SM_${NORMAL}.${CHR}.DR.sort.bam -normal ${NORMAL}.${CHR}.DR.bam -O ${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz 
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz -O filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_name.vcf.gz
    # sample name and PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SM_${TUMOR}.${CHR}.DR.sort.bam -tumor ${TUMOR}.${CHR}.DR.bam -I SM_${NORMAL}.${CHR}.DR.sort.bam -normal ${NORMAL}.${CHR}.DR.bam -O ${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz --panel-of-normals  ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz -O filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_name_PON.vcf.gz
    # PON
    gatk Mutect2 -R ${PATH_REF}${REF} -I SM_${TUMOR}.${CHR}.DR.sort.bam -I SM_${NORMAL}.${CHR}.DR.sort.bam -O ${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz --panel-of-normals  ${PON}
    gatk FilterMutectCalls -R ${PATH_REF}${REF} -V ${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz -O filtered_${TUMOR}_${NORMAL}_${CHR}_somatic_PON.vcf.gz
}

Filter_with_hand() {

}