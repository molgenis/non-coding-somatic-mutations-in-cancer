'''
Created on 12 feb. 2022

@author: roy.oelen
'''

import operator
import numpy

class VcfMatrix:
    '''
    classdocs
    '''


    def __init__(self, participants=None, snps=None, matrix=None, donor_name_mapping=None):
        '''
        Constructor
        '''
        self.participants = participants
        self.snps = snps
        self.matrix = matrix
        self.donor_name_mapping = donor_name_mapping

    # override to be more informative
    def __str__(self):
        '''
        get a string representation of the object
        :return: a string representation of the object
        '''
        #message = ''.join(['genotype matrix across ' , str(self.participants.len()), ' participants and ' , str(self.snps.len()) , ' snps\n'])
        message = 'banaan'
        return message

    def get_snps(self):
        return self.snps

    def get_participants(self):
        return self.participants

    def get_matrix(self):
        return self.matrix

    def set_order(self, order):
        # create a new matrix
        matrix_copy = numpy.empty((len(order), self.participants.size), dtype='float64')
        # empty the matrix
        matrix_copy[:] = numpy.NaN
        # start filling
        for index_old in range(0, len(order), 1):
            # get the new index
            index_new = order[index_old]
            matrix_copy[index_old] = self.matrix[index_new]
        # replace the old matrix
        self.matrix = None
        self.matrix = matrix_copy
        # set the SNP order
        self.snps = self.snps[order]

    def set_donor_name_mapping(self, donor_name_mapping):
        self.donor_name_mapping = donor_name_mapping


    def get_donor_names(self, donor_ids):
        return operator.itemgetter(donor_ids)(self.donor_name_mapping)


    def get_donor_names(self, donor_ids):
        # initialize
        donor_names = []
        # check each donor id
        for donor_id in donor_ids:
            # add the correct donor name
            donor_name = self.donor_name_mapping[donor_id]
            # add to the list
            donor_names.append(donor_name)
        return donor_names