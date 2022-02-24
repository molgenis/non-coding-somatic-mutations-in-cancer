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
codeList=()


ml BCFtools/1.11-GCCcore-7.3.0
 

for filename in ${PATH_OLD}*.vcf.gz*; do
    # Get basename
    BASENAME=$( echo "$(basename -- $filename)")
    # Get code (before .)
    IFS='.' read -r codeName string <<< "$BASENAME"
    # Check if codeName (=code) exist in codeList
    # If it is not in the list it is added and the files are combined with that code.
    if [[ ! " ${codeList[*]} " =~ " ${codeName} " ]]; then
        # Index files
        bcftools index ${PATH_OLD}${codeName}.wgs.CaVEMan.raw_somatic_mutation.vcf.gz
        bcftools index ${PATH_OLD}${codeName}.wgs.sanger_raw_pindel.raw_somatic_mutation.vcf.gz
        # Concat files
        bcftools concat -a -Oz ${PATH_OLD}${codeName}.wgs.CaVEMan.raw_somatic_mutation.vcf.gz ${PATH_OLD}${codeName}.wgs.sanger_raw_pindel.raw_somatic_mutation.vcf.gz > ${PATH_NEW}${sepName}.vcf.gz
        # echo 'JA'
        # echo ${codeName}
        codeList+=( ${codeName} )
    fi
    
done

# echo ${codeList[*]}