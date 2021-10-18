 for f in *.DR.bam; do samtools view -H $f  | sed "s/SM:[^\t]*/SM:${f}/g" | samtools reheader - $f > SM_$f; echo -ne "$f\t" ; samtools view -H $f | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq; done


for f in *.DR.bam; do echo -ne "$f\t" ; samtools view -H $f  | sed "s/SM:[^\t]*/SM:${f}/g" | samtools reheader - $f > SM_$f; done





for f in *.DR.bam; do echo -ne "$f\t" ; samtools view -H $f  | sed "s/SM:[^\t]*/SM:${f}/g" | samtools reheader - $f > SM_$f; done
for f in *.DR.bam; do echo -ne "$f\t" ; samtools view -H $f | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq; done


SM_SS6005042.chr22.DR.bam       SS6005042.chr22.DR.bam
SM_SS6005043.chr22.DR.bam       SS6005043.chr22.DR.bam
SM_SS6005044.chr22.DR.bam       SS6005044.chr22.DR.bam
SS6005042.chr22.DR.bam  sam_id
SS6005043.chr22.DR.bam  sam_id
SS6005044.chr22.DR.bam  sam_id


gatk Mutect2 -R /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa -I SM_SS6005044.chr22.DR.bam -tumor SS6005044.chr22.DR.bam -I SM_SS6005042.chr22.DR.bam -normal SS6005042.chr22.DR.bam -O 5044_5042_somatic_name.vcf.gz


gatk Mutect2 -R /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa -I SM_SS6005042.chr22.DR.bam -normal SS6005042.chr22.DR.bam -O 5042_somatic_name.vcf.gz


gatk Mutect2 -R /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/chr/unzip/chr22.fa -I SM_SS6005044.chr22.DR.sort.bam -I SM_SS6005043.chr22.DR.sort.bam -I SM_SS6005042.chr22.DR.sort.bam -normal SS6005042.chr22.DR.bam -O 5044_5043_5042_somatic_name.vcf.gz

--panel-of-normals /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/PanelOfNormals/somatic-b37_Mutect2-WGS-panel-b37.vcf