'''
Created on 12 feb. 2022

@author: roy.oelen
'''

import VcfInfo
import VcfMatrix
import numpy

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

    def harmonize(self):
        '''
        harmonize the object as such that the VcfInfo and the VcfMatrix only contain snps in both objects, and that the order is the same in both
        :return:
        '''
        # get the SNPs we definitely have data for
        snps_info = self.vcf_info.snp_ids
        snps_geno = self.vcf_matrix.snps
        snps_info = set(snps_info)
        snps_geno = set(snps_geno)
        snps_common = snps_info.intersection(snps_geno)
        # sort these
        snps_common = sorted(snps_common)
        # initialize positions of the SNPs
        info_locs = numpy.empty((len(snps_common)), dtype='intc')
        matrix_locs = numpy.empty((len(snps_common)), dtype='intc')
        # empty the matrix
        #info_locs[:] = numpy.NaN
        #matrix_locs[:] = numpy.NaN
        info_locs[:] = -1
        matrix_locs[:] = -1
        # now do this from beginning to end
        for i in range(0,len(snps_common), 1):
            # this is the SNP
            snp = snps_common[i]
            # this is the location in the info
            info_loc = numpy.argmax(self.vcf_info.snp_ids == snp)
            # this is the genotype loc
            geno_loc = numpy.argmax(self.vcf_matrix.snps == snp)
            # set these locations in the location array
            info_locs[i] = info_loc
            matrix_locs[i] = geno_loc
        # now set the order in the info
        self.vcf_info.set_order(info_locs)
        # and set that order in the genotype matrix
        self.vcf_matrix.set_order(geno_loc)
        print(''.join(['harmonised ', str(self)]))

    def get_donor_names(self, donor_ids):
        return self.get_donor_names(donor_ids)