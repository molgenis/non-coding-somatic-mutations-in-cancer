'''
Created on 12 feb. 2022

@author: roy.oelen
'''

import deprecation
import pandas

class VcfInfo:
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.snps = None
        self.info = None


    def set_data_columns(self, chromosomes, positions, identifiers, references, alts, qualts, filters, infos):
        '''
        set the info data for the object

        :param chromosomes: a numpy array of chromosomes
        :param positions: a numpy array of the positions on the chromosome
        :param identifiers: an identier for each entry in a numpy array
        :param references: a reference allele for each entry in a numpy array
        :param alts: a alternative allele for each entry in a numpy array
        :param qualts: a numpy array that is not used yet
        :param filters: a numpy array that is not used yet
        :param infos: a numpy array that is not used yet
        :return: NONE the object is modified in-place
        '''

        # create the dataframe
#        data = pandas.DataFrame({
#            '#CHROM' : chromosomes,
#            'POS' : positions,
#            'ID' : identifiers,
#            'REF' : references,
#            'ALT' : alts,
#            'QUAL' : qualts,
#            'FILTER' : filters,
#            'INFO' : infos})

        data = {
            '#CHROM' : chromosomes,
            'POS' : positions,
            'ID' : identifiers,
            'REF' : references,
            'ALT' : alts,
            'QUAL' : qualts,
            'FILTER' : filters,
            'INFO' : infos}
        self.info = data

    @deprecation.deprecated(details="now using numpy for efficiency sake")
    def set_data_pandas(self, info_panda):
        '''

        '''
        # copy a dataframe
        data = info_panda
        self.info = info_panda
        self.snps = info_panda.index.values()