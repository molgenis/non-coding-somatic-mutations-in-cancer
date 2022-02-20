'''
Created on 12 feb. 2022

@author: roy.oelen
'''

# external imports
import pandas
import numpy
# internal imports
import VcfInfo
import VcfMatrix

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
    def tuples_to_vcfdata(list_of_tuples, snp_ids, donor_ids, dosage=True, lower=0.1, higher=0.9):
        # get the number of snps and ids
        nr_snp_ids = len(snp_ids)
        nr_donor_ids = len(donor_ids)
        # we'll fill by coordinates in the 2d array, so we map ids to indices
        snp_zip = zip(snp_ids, range(0, nr_snp_ids+1, 1))
        snp_id_coordinates = dict(snp_zip) # gives you this {'snp_id1' : 0, 'snp_id2' : 1, 'snp_id3 : 2'}
        donor_zip = zip(donor_ids, range(0, nr_donor_ids+1, 1))
        donor_id_coordinates = dict(donor_zip)
        # reserve the memory
        gt_matrix = numpy.empty((nr_snp_ids, nr_donor_ids))
        # empty the matrix
        gt_matrix[:] = numpy.NaN
        # we'll use this opportunity to also determine the order of the donors and genes
        sorted_donor_ids = numpy.empty((1, nr_donor_ids))
        sorted_snp_ids = numpy.empty((1, nr_snp_ids))

        # query_start = 'SELECT dhs.snp_ID, dhs.donor_ID, dhs.mutant_allele_read_count, dhs.total_read_count
        for tuple_entry in list_of_tuples:
            print(tuple_entry)
            # get the snp ID entry in the tuple
            #snp_id = tuple_entry['snp_ID']
            snp_id = tuple_entry[0]
            # get the coordinate of the entry
            snp_pos = snp_id_coordinates[snp_id]
            # same for donor
            #donor_id = tuple_entry['donor_ID']
            donor_id = tuple_entry[1]
            donor_pos = donor_id_coordinates[donor_id]
            # get the counts
            #alt_reads = tuple_entry['mutant_allele_read_count']
            alt_reads = tuple_entry[2]
            #all_reads = tuple_entry['total_read_count']
            all_reads = tuple_entry[3]
            # init variable
            dosage = None
            # check the possible reads
            if all_reads is None:
                dosage = None
            elif all_reads is not None and alt_reads is None:
                dosage = 0
            elif all_reads is not None and alt_reads is not None:
                # calculate the dosage
                dosage = alt_reads / all_reads

            # either do dosage or convert to genotypes
            if dosage is False:
                gt_matrix[snp_pos][donor_pos] = dosage
            # taking into account that dosage might be empty
            elif dosage is True and dosage is not None:
                if dosage < lower:
                    gt_matrix[snp_pos][donor_pos] = 0
                elif dosage >= lower and dosage < higher:
                    gt_matrix[snp_pos][donor_pos] = 1
                elif dosage >= higher:
                    gt_matrix[snp_pos][donor_pos] = 2
            # we know that in this location, there IDs are present, so let's save them
            sorted_snp_ids[0][snp_pos] = snp_id
            sorted_donor_ids[0][donor_pos] = donor_id
        # now turn it into a VcfMatrix object
        vcf_matrix = VcfMatrix(sorted_donor_ids, sorted_snp_ids, gt_matrix)

        return vcf_matrix