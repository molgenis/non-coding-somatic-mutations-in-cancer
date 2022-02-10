#!/usr/bin/env python3
from SNP import SNP

class OVERVIEW_SNP:

    def __init__(self):
        """
        Constructor
        """
        self.vcf_dict = dict()
        self.dict_SNP_ID = dict()
        self.all_donors = dict()

        

    def set_snp(self, chr, pos, ref, alt, donor_id, total_read_count, 
                mutant_allele_read_count, icgc_specimen_id, project_code):
        
        SNP_ID = f'{chr}_{pos}_{ref}_{alt}'
        homo_hetero, total_read_count, mutant_allele_read_count = self.check_homo_hetero(total_read_count, mutant_allele_read_count)
        snp_format_donor = ':'.join([str(homo_hetero), str(total_read_count), str(mutant_allele_read_count), str(icgc_specimen_id), str(project_code)]) #[icgc_specimen_id, project_code, homo_hetero, total_read_count, mutant_allele_read_count]
        
        if donor_id not in self.all_donors:
            self.all_donors[donor_id] = ':'.join([str('0/0'), str('.'), str('.'), str(icgc_specimen_id), str(project_code)])
        
        if SNP_ID not in self.dict_SNP_ID:
            snp = SNP(chr, pos, ref, alt)
            self.dict_SNP_ID[SNP_ID] = snp
        else:
            snp = self.dict_SNP_ID[SNP_ID]
        

        if SNP_ID in self.vcf_dict:
            if donor_id not in self.vcf_dict[SNP_ID]:
                self.vcf_dict[SNP_ID][donor_id] = snp_format_donor
        else:
            # new_donor_dict = dict()
            # new_donor_dict[snp.donor_id] = [icgc_specimen_id, project_code, snp.homo_hetero]
            # self.vcf_dict[SNP_ID] = {'SNP_INFO': snp, donor_id: snp_format_donor}
            self.vcf_dict[SNP_ID] = {'CHROM': snp.chr, 'POS': snp.pos, 'ID': snp.id, 'REF': snp.ref, 'ALT': snp.alt,
                                    'QUAL': snp.qual, 'FILTER': snp.filter, 'INFO': snp.info, 'FORMAT': snp.format,
                                    donor_id: snp_format_donor}

    def check_homo_hetero(self, total_read_count, mutant_allele_read_count):
        if mutant_allele_read_count != '' and total_read_count != '':
            if int(mutant_allele_read_count)/int(total_read_count) >= 0.9:
                homo_hetero = '1/1'
            else:
                homo_hetero = '0/1'
        else:
            if mutant_allele_read_count == '':
                mutant_allele_read_count = '.'
            if total_read_count == '':
                total_read_count = '.'
            homo_hetero = '.'
        return homo_hetero, total_read_count, mutant_allele_read_count

    def get_vcf_dict(self):
        return self.vcf_dict

    def get_all_donors(self):
        return self.all_donors


    
        


    
        

    







