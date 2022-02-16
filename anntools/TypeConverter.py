'''
Created on 12 feb. 2022

@author: roy.oelen
'''

import pandas
import numpy
import VcfInfo

class TypeConverter:


    @staticmethod
    def sqlite3_to_pandas(sqlite3_output, columns = None):
        '''
        convert tuple from sqlite3 output to a pandas dataframe
        :param sqlite3_output: tuples from sqlite3 output
        :param columns: what to name columns
        :return: a pandas dataframe
        '''
        pandas_df = pandas.DataFrame(sqlite3_output, columns)
        return pandas_df

    @staticmethod
    def pandas_to_numpy(pandas_df):
        '''
        convert pandas dataframe to numpy array
        :param pandas_df: the pandas dataframe to convert
        :return: the resulting numpy array
        '''
        numpy_array = pandas_df.to_numpy()
        return numpy_array

    @staticmethod
    def pandas_to_vcfcompatible(pandas_df):
        # TODO actually implement this
        print('pandas_to_vcfcompatible not yet implemented')
        return None

    @staticmethod
    def tuples_to_vcfinfo(list_of_tuples):
        '''
        convert tuples like from sqlite3 output to a VcfInfo object
        :param list_of_tuples: the list of tuples to convert
        :return:  a VcfInfo object
        '''
        # TODO make this not depend on the order
        snp_IDs, snp_chrs, snp_pos_starts, snp_pos_ends, snp_refs, snp_alts = zip(*list_of_tuples)
        # make identifiers and other columns that are currently empty
        identifiers = numpy.empty([1, len(list_of_tuples)],dtype='U')
        qualts = numpy.empty([1, len(list_of_tuples)])
        filters = numpy.empty([1, len(list_of_tuples)])
        infos = numpy.empty([1, len(list_of_tuples)])
        # fill everything
        for i in range(0, len(list_of_tuples), 1):
            qualts[0][i] = None
            filters[0][i] = None
            infos[0][i] = None
            qualts[0][i] = None
            snp_chr = str(snp_chrs[i])
            snp_pos = str(snp_pos_starts[i])
            identifiers[0][i] = ':'.join([snp_chr, snp_pos])

        # convert to numpy where needed
        snp_chrs = numpy.array(snp_chrs)
        snp_pos_starts = numpy.array(snp_pos_starts)
        # create a new VcfInfo
        vcf_info = VcfInfo.VcfInfo()
        vcf_info.set_data_columns(snp_chrs,snp_pos_starts,identifiers, snp_refs,snp_alts,qualts,filters,infos)
        return vcf_info


    @staticmethod
    def tuples_to_vcfdata(list_of_genos):
        # TODO actually implement this
        print('tuples_to_vcfdata not yet implemented')
        return None