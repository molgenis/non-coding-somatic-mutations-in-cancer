'''
Created on 12 feb. 2022

@author: roy.oelen
'''

import operator


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
        # set the new order of the matrix
        self.matrix = self.matrix[order]
        # set the SNP order
        self.snps = self.snps[order]

    def set_donor_name_mapping(self, donor_name_mapping):
        self.donor_name_mapping = donor_name_mapping


    def get_donor_names(self, donor_ids):
        return operator.itemgetter(donor_ids)(self.donor_name_mapping)