#!/usr/bin/env python3

class SNP:

    def __init__(self, chr, pos, ref, alt):
        """
        Constructor
        """
        self.chr = chr
        self.pos = pos
        self.id = '.'
        self.ref = ref
        self.alt = alt
        self.qual = '.'
        self.filter = '.'
        self.info = '.'
        # self.SNP_ID = f'{chr}_{pos}_{ref}_{alt}'
        # self.donor_id = donor_id
        # self.total_read_count = total_read_count
        # self.mutant_allele_read_count = mutant_allele_read_count
        # self.homo_hetero = self.check_homo_hetero(total_read_count, mutant_allele_read_count)
        # self.icgc_specimen_id = icgc_specimen_id
        # self.project_code = project_code
        
        # GT = homo_hetero, TC = total_read_count, MC = mutant_allele_read_count, SP = icgc_specimen_id, PRC = project_code
        self.format = 'GT:TC:MC:SP:PRC'


    
        
        

    







