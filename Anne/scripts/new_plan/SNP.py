#!/usr/bin/env python3

class SNP:

    def __init__(self, chr, pos, ref, alt, project_code):
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
        self.info = f'project_code={project_code}'      
        # GT = homo_hetero, TC = total_read_count, MC = mutant_allele_read_count,
        # SP = specimen_id, PRC = project_code
        self.format = 'GT:TC:MC:SP'


    
        
        

    







