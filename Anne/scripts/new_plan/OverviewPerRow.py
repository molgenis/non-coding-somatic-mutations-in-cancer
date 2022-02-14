#!/usr/bin/env python3
from SNP import SNP


class OverviewPerRow:

    def __init__(self):
        """
        Constructor
        """
        self.vcf_dict = dict()
        self.dict_SNP_ID = dict()
        self.all_donors = dict({'CHROM':'', 'POS':'', 'ID':'', 'REF':'', 'ALT':'', 'QUAL':'', 'FILTER':'', 'INFO':'', 'FORMAT':'',})

    def set_snp(self, chr, pos, ref, alt, donor_id, total_read_count,
                mutant_allele_read_count, specimen_id, project_code):
        """

        """
        # Create snp ID
        SNP_ID = f'{chr}_{pos}_{ref}_{alt}'
        # Call check_homo_hetero
        homo_hetero, total_read_count, mutant_allele_read_count = self.check_homo_hetero(total_read_count,
                                                                                         mutant_allele_read_count, project_code) #TODO DELETE
        # Create format structure with values
        snp_format_donor = ':'.join(
            [str(homo_hetero), str(total_read_count), str(mutant_allele_read_count), str(specimen_id)])
        # [specimen_id, project_code, homo_hetero, total_read_count, mutant_allele_read_count]
        # Make a dictionary with the donor ids as keys and the format structure as values,
        # but then for the donors who do not have that snp.
        # So they don't have homo_hetero, total_read_count, mutant_allele_read_count.
        if donor_id not in self.all_donors:
            self.all_donors[donor_id] = ':'.join(
                [str('.'), str('.'), str('.'), str(specimen_id)])

        # Check if the snp id is in the vcf_dict. 
        if SNP_ID in self.vcf_dict:
            # Check whether the donor_id already exists in vcf_dict[SNP_ID], 
            # if not add this donor (with the format structure of this donor).
            if donor_id not in self.vcf_dict[SNP_ID]:
                self.vcf_dict[SNP_ID][donor_id] = snp_format_donor
        # If it doesn't exist, create the snp object.
        # And add value to vcf_dict[SNP_ID]. vcf_dict[SNP_ID] is a dictionary that consists of the 
        # standard columns of a vcf file as keys but also the donors with their format structure.
        else:
            snp = SNP(chr, pos, ref, alt, project_code)
            self.vcf_dict[SNP_ID] = {'CHROM': snp.chr, 'POS': snp.pos, 'ID': snp.id, 'REF': snp.ref, 'ALT': snp.alt,
                                     'QUAL': snp.qual, 'FILTER': snp.filter, 'INFO': snp.info, 'FORMAT': snp.format,
                                     donor_id: snp_format_donor}

    def check_homo_hetero(self, total_read_count, mutant_allele_read_count, project_code):
        """

        """
        # Check if mutant_allele_read_count and total_read_count are empty
        # If these are not empty divide mutant_allele_read_count by total_read_count
        # and determine if they are homozygous or heterozygous.
        if mutant_allele_read_count != '' and total_read_count != '':
            if int(mutant_allele_read_count) != 0 or int(total_read_count) != 0:
                if int(mutant_allele_read_count) / int(total_read_count) >= 0.9:
                    homo_hetero = '1/1'
                else:
                    homo_hetero = '0/1'
            else:
                print(f'mutant_allele_read_count {mutant_allele_read_count}')
                print(f'total_read_count {total_read_count}')
                print(f'project_code {project_code}')
                print('-----')
                homo_hetero = '.'

        # If mutant_allele_read_count and total_read_count are empty,
        # change both variables with a '.'.
        # Make homo_hetero a '.' for now. # TODO TODO
        else:
            if mutant_allele_read_count == '':
                mutant_allele_read_count = '.'
            if total_read_count == '':
                total_read_count = '.'
            homo_hetero = '.'
        return homo_hetero, total_read_count, mutant_allele_read_count

    def get_vcf_dict(self):
        """

        """
        return self.vcf_dict

    def get_all_donors(self):
        """

        """
        return self.all_donors
