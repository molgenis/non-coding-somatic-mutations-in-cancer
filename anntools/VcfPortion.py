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


    def __init__(self, vcf_matrix=None, vcf_info=None):
        '''
        Constructor
        '''
        self.vcf_matrix = vcf_matrix
        self.vcf_info = vcf_info


    def __str__(self):
        message = ''.join([ 'VcfPortion object\n' , 'VcfInfo\n' , str(self.vcf_info) , '\n' , 'VcfFmatrix\n' + str(self.vcf_matrix) , '\n'])
        return message

    def get_snps_by_ids(self, identifiers):
        print('get_snps_by_ids not yet implemented')
        return None
