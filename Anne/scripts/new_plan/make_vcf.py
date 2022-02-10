import gzip
from OVERVIEW_SNP import OVERVIEW_SNP
from SNP import SNP

project_file = gzip.open("D:/Hanze_Groningen/STAGE/NEW PLAN/ALL-US.tsv.gz",'rt')
header = project_file.readline().strip().split("\t")

overview = OVERVIEW_SNP()

for line in project_file:
    line = line.strip()
    elems = line.split("\t")
    if elems[33] == 'WGS':
        donor_id = elems[1]
        project_code = elems[2]
        icgc_specimen_id = elems[3]
        chr = elems[8]
        pos = elems[9]
        ref = elems[15]
        alt = elems[16]
        total_read_count = elems[19]
        mutant_allele_read_count = elems[20]
        overview.set_snp(chr, pos, ref, alt, donor_id, total_read_count, 
                        mutant_allele_read_count, icgc_specimen_id, project_code)


print(overview.get_vcf_dict())
print(len(overview.get_all_donors()))

vcf_dict = overview.get_vcf_dict()
all_donors = overview.get_all_donors()

for SNP_ID, vcf_info in vcf_dict.items():
    for donor_ID, format_info in all_donors.items():
        if donor_ID not in vcf_info:
            vcf_dict[SNP_ID][donor_ID] = format_info

for key, value in vcf_dict.items():
    if len(value) != (len(all_donors) + 9):
        print(f'--------------{key}')

print(overview.get_vcf_dict())





