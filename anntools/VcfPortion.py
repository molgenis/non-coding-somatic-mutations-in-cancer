'''
Created on 12 feb. 2022

@author: roy.oelen
'''

import VcfInfo
import VcfMatrix
import pandas


class VcfPortion(object):
    '''
    classdocs
    '''


    def __init__(self, genotype_matrix, vcf_info):
        '''
        Constructor
        '''
        self.genotype_matrix = None
        self.vcf_info = None


    def get_snps_by_ids(self, identifiers):
        # subset the info portion
        subset_info_portion = self.vcf_info[self.vcf_info['ID'].isin(identifiers)]
        # create a new info portion
        subset_vcf_info = VcfInfo()
        subset_vcf_info.set_data_panda(subset_info_portion)

        # subset the matrix portion
        subset_matrix_portion = self.vcf_info[self.vcf_info['ID'].isin(identifiers)]
        # create a new matrix
        subset_vcf_matrix = VcfMatrix()
        subset_vcf_matrix.set_data_pandas(subset_matrix_portion)

        # create the object
        vcfPortion = VcfPortion(subset_vcf_matrix, subset_vcf_info)
        # TODO finish 
        print('get_snps_by_ids not yet implemented')
