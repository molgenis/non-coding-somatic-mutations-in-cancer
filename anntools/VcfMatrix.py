'''
Created on 12 feb. 2022

@author: roy.oelen
'''

class VcfMatrix:
    '''
    classdocs
    '''


    def __init__(self, participants=None, snps=None, matrix=None):
        '''
        Constructor
        '''
        self.participants = participants
        self.snps = snps
        self.matrix = matrix

    # override to be more informative
    def __str__(self):
        '''
        get a string representation of the object
        :return: a string representation of the object
        '''
        message = ''.join(['genotype matrix across ' , str(len(self.participants)), ' participants and ' , str(len(self.snps)) + ' snps\n'])
        return message


