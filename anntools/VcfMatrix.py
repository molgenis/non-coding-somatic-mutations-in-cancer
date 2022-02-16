'''
Created on 12 feb. 2022

@author: roy.oelen
'''

class VcfMatrix:
    '''
    classdocs
    '''


    def __init__(self, params):
        '''
        Constructor
        '''
        self.participants = None
        self.snps = None
        self.matrix = None


    def set_data_pandas(self, matrix_panda):
        '''

        '''
        # copy a dataframe
        self.matrix = matrix_panda
        self.snps = matrix_panda.index.values()