'''
Created on 12 feb. 2022

@author: roy.oelen
'''

import VcfInfo
import VcfMatrix
import numpy
import math

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
        self.is_sorted = False


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
        self.vcf_matrix.set_order(matrix_locs)
        print(''.join(['harmonised ', str(self)]))
        # set the hormonized tag
        self.is_sorted = True

    def get_donors(self, names=False):
        # fetch the participants
        donor_ids = self.vcf_matrix.get_participants()
        # return them if we don't need to get the names
        if names is False or names is None:
            return(donor_ids)
        else:
            donor_names = self.get_donor_names(donor_ids)
            return donor_names

    def get_donor_names(self, donor_ids):
        return self.vcf_matrix.get_donor_names(donor_ids)

    def subset_self(self, nr_horizontal_slices=2, nr_vertical_slices=2):
        # get the number of snps, the horizontal slice
        nr_snps = self.vcf_info.snp_ids.size
        # get the number of donors, the vertical slice
        nr_donors =  self.vcf_matrix.participants.size
        # calculate the total number of slices
        total_nr_slices = nr_horizontal_slices * nr_vertical_slices
        # reserve that amount of space
        slices = [None] * total_nr_slices
        # calculate the size of the slizes
        snp_slice_size = math.ceil(nr_snps / nr_horizontal_slices)
        donor_slice_size = math.ceil(nr_donors / nr_vertical_slices)


    def get_slice(self, vertical_start, vertical_stop, horizontal_start, horizontal_stop):
        # we need to make sure that the positions are correct
        if self.is_sorted is False:
            self.harmonize()
        # grab the donors
        donor_slice_ids = self.vcf_matrix.get_participants()[horizontal_start : horizontal_stop]
        # grab the SNPs
        snp_slice_ids = self.vcf_matrix.snps[vertical_start : vertical_stop]
        # grab the matrix
        matrix_slice = self.vcf_matrix.matrix[vertical_start : vertical_stop][horizontal_start : horizontal_stop]
