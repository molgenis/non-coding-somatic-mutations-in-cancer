'''
Created on 12 feb. 2022

@author: roy.oelen
'''

import numpy
class VcfInfo:
    '''
    classdocs
    '''


    def __init__(self, snp_ids=None, chromosomes=None, positions=None, identifiers=None, references=None, alts=None, qualts=None, filters=None, infos=None):
        '''
        #TODO documentations here
        :param chromosomes:
        :param positions:
        :param identifiers:
        :param references:
        :param alts:
        :param qualts:
        :param filters:
        :param infos:
        '''
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
        self.snp_ids = snp_ids

    # override to be more informative
    def __str__(self):
        '''
        get a string representation of the object
        :return: a string representation of the object
        '''
        message = ''.join(['instance of VcfInfo across ' , str(len(self.info['ID'][0])) , ' genotypes, at ' , str(len(numpy.unique(self.info['#CHROM']))) , ' chromosomes, and ' , str(len(numpy.unique(self.info['POS']))) , ' positions\n'])
        print(str(self.info['ID']))
        return message

    def set_data_columns(self, snp_ids, chromosomes, positions, identifiers, references, alts, qualts, filters, infos):
        '''
        set the info data for the object

        :param snp_ids: a numpy array of IDs as they are in the database
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
        self.snp_ids = snp_ids

    def get_info(self):
        return self.info


    def set_order(self, order):
        '''

        :param order:
        :return:
        '''
        # get the keys
        keys = self.info.keys()
        # go through all of the keys data, and set the order
        for key in keys:
            # grab data
            value = self.info[key]
            # check if the value is anything
            if value is not None:
                # order this numpy array
                new_value = numpy.take(value, order)
                # set the new value
                # TODO check if this is safe!
                self.info[key] = new_value
        # set the snp order as well
        self.snp_ids = self.snp_ids[order]

