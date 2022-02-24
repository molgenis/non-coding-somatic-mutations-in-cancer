#!/usr/bin/bash

#SBATCH --job-name=concat
#SBATCH --output=concat.out
#SBATCH --error=concat.err
#SBATCH --time=49:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

PATH_NEW=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/merge_D/
PATH_OLD=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/All_B/
# PATH_OLD=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GDC/All_H/
listje=()


ml BCFtools/1.11-GCCcore-7.3.0
 

for filename in ${PATH_OLD}*.vcf.gz*; do
    #echo "$filename"
    BASENAME=$( echo "$(basename -- $filename)")
    #echo ${BASENAME}
    #sepName=$(echo $BASENAME | tr ".wgs" "\n")
    IFS='.' read -r sepName string <<< "$BASENAME"
    if [[ ! " ${listje[*]} " =~ " ${sepName} " ]]; then
        bcftools index ${PATH_OLD}${sepName}.wgs.CaVEMan.raw_somatic_mutation.vcf.gz
        bcftools index ${PATH_OLD}${sepName}.wgs.sanger_raw_pindel.raw_somatic_mutation.vcf.gz
        bcftools concat -a -Oz ${PATH_OLD}${sepName}.wgs.CaVEMan.raw_somatic_mutation.vcf.gz ${PATH_OLD}${sepName}.wgs.sanger_raw_pindel.raw_somatic_mutation.vcf.gz > ${PATH_NEW}${sepName}.vcf.gz
        echo 'JA'
        echo ${sepName}
        listje+=( ${sepName} )
    fi
    
done

echo ${listje[*]}