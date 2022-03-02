'''
Created on 12 feb. 2022

@author: roy.oelen
'''

# external imports
import pandas
import numpy
from io import StringIO
# internal imports
import VcfInfo
import VcfMatrix
import StringBuilder

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
        identifiers = numpy.empty([len(list_of_tuples)],dtype='<U12')
        qualts = numpy.empty([len(list_of_tuples)])
        filters = numpy.empty([len(list_of_tuples)])
        infos = numpy.empty([len(list_of_tuples)])
        # dirty workaround
        #ids = list()
        #snp_refs = numpy.empty([1, len(list_of_tuples)])
        #snp_alts = numpy.empty([1, len(list_of_tuples)])
        # fill everything
        for i in range(0, len(list_of_tuples), 1):
            qualts[i] = None
            filters[i] = None
            infos[i] = None
            qualts[i] = None
            snp_chr = str(snp_chrs[i])
            snp_pos = str(snp_pos_starts[i])
            ident = ':'.join([snp_chr, snp_pos])
            identifiers[i] = ident
            #ids.append(ident)

        # convert to numpy where needed
        snp_chrs = numpy.array(snp_chrs)
        snp_pos_starts = numpy.array(snp_pos_starts, dtype='intc')
        # ref/alt as well
        snp_refs = numpy.array(snp_refs)
        snp_alts = numpy.array(snp_alts)
        # and the IDs like in the database
        snp_ids = numpy.array(snp_IDs, dtype='intc')
        # create a new VcfInfo
        #vcf_info = VcfInfo.VcfInfo(snp_chrs, snp_pos_starts, numpy.array(ids), snp_refs, snp_alts, None, None, None)
        vcf_info = VcfInfo.VcfInfo(snp_ids, snp_chrs, snp_pos_starts, identifiers, snp_refs, snp_alts, None, None, None)
        return vcf_info


    @staticmethod
    def tuples_to_vcfdata(list_of_tuples, snp_ids, donor_ids, use_dosage=True, lower=0.1, higher=0.9):
        # get the number of snps and ids
        nr_snp_ids = len(snp_ids)
        nr_donor_ids = len(donor_ids)
        # we'll fill by coordinates in the 2d array, so we map ids to indices
        snp_zip = zip(snp_ids, range(0, nr_snp_ids+1, 1))
        snp_id_coordinates = dict(snp_zip) # gives you this {'snp_id1' : 0, 'snp_id2' : 1, 'snp_id3 : 2'}
        donor_zip = zip(donor_ids, range(0, nr_donor_ids+1, 1))
        donor_id_coordinates = dict(donor_zip)
        # reserve the memory
        gt_matrix = numpy.empty((nr_snp_ids, nr_donor_ids), dtype='float64')
        # empty the matrix
        gt_matrix[:] = numpy.NaN

        # we'll use this opportunity to also determine the order of the donors and genes
        sorted_donor_ids = numpy.empty(nr_donor_ids, dtype='intc')
        sorted_snp_ids = numpy.empty(nr_snp_ids, dtype='intc')
        #sorted_donor_ids[:] = numpy.NaN
        #sorted_snp_ids[:] = numpy.NaN
        sorted_donor_ids[:] = -1
        sorted_snp_ids[:] = -1

        # query_start = 'SELECT dhs.snp_ID, dhs.donor_ID, dhs.mutant_allele_read_count, dhs.total_read_count
        for tuple_entry in list_of_tuples:
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
            else:
                print(' '.join(['bad entry for', donor_id, snp_id]))
            # either do dosage or convert to genotypes
            if use_dosage is True:
                gt_matrix[snp_pos][donor_pos] = dosage
            # taking into account that dosage might be empty
            elif use_dosage is False and dosage is not None:
                if dosage < lower:
                    gt_matrix[snp_pos][donor_pos] = 0
                elif dosage >= lower and dosage < higher:
                    gt_matrix[snp_pos][donor_pos] = 1
                elif dosage >= higher:
                    gt_matrix[snp_pos][donor_pos] = 2
            else:
                print(' '.join(['bad entry for', str(donor_id), str(snp_id)]))
            # we know that in this location, there IDs are present, so let's save them
            sorted_snp_ids[snp_pos] = int(snp_id)
            sorted_donor_ids[donor_pos] = int(donor_id)
        # now turn it into a VcfMatrix object
        vcf_matrix = VcfMatrix.VcfMatrix(sorted_donor_ids, sorted_snp_ids, gt_matrix)

        return vcf_matrix

    @staticmethod
    def vcf_portion_to_file(vcf_portion, full_output_loc):
        # get the donors
        donor_names = vcf_portion.get_donors(names=True)
        # get the SNPs
        TypeConverter.write_vcf_content(vcf_portion, full_output_loc, donors = donor_names)


    @staticmethod
    def make_vcf_header(vcf_portion, donor_name=True, donors=None):
        # the first part is always the same
        header_start = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
        # get the donors
        donor_names = donors
        # if they were not supplied get them yourself
        if donor_names is None:
            donor_names = vcf_portion.get_donors(names=True)

        # paste them together
        donor_line = '\t'.join(donor_names)
        # paste together the complete line
        header = '\t'.join([header_start, donor_line])
        return header

    @staticmethod
    def write_vcf_content(vcf_portion, full_output_loc, header=None, donors=None):
        # first create the header
        header_to_use = header
        if header_to_use is None:
            # make a header if none is supplied
            header_to_use = TypeConverter.make_vcf_header(vcf_portion, donor_name=True)

        # start to open a file
        file = open(full_output_loc, 'a')

        # we will loop through the donor_ids and snp_ids
        #for snp in numpy.nditer(vcf_portion.vcf_matrix, flags = ['refs_ok']):
        for index_snp in range(0, vcf_portion.vcf_matrix.snps.size, 1):
            # extract the snp name
            snp_name = vcf_portion.vcf_info.snp_ids[index_snp]

            # get the metadata of the SNP
            chrom = vcf_portion.vcf_info.get_info('#CHROM')[index_snp]
            pos = vcf_portion.vcf_info.get_info('POS')[index_snp]
            id = vcf_portion.vcf_info.get_info('ID')[index_snp]
            ref = vcf_portion.vcf_info.get_info('REF')[index_snp]
            alt = vcf_portion.vcf_info.get_info('ALT')[index_snp]
            qual = vcf_portion.vcf_info.get_info('QUAL')[index_snp]
            filter = vcf_portion.vcf_info.get_info('FILTER')[index_snp]
            info = vcf_portion.vcf_info.get_info('INFO')[index_snp]
            # append them to the line in the vcf file
            string_builder = StringBuilder.StringBuilder()
            string_builder.append('\t'.join([chrom, str(pos), str(snp_name), ref, alt, str(qual), str(filter), str(info)]))

            # check each donor
            #for call in numpy.nditer(snp, flags = ['refs_ok']):
            for index_donor in range(0, len(vcf_portion.vcf_matrix.participants), 1):
                # get the call
                #string_builder.append('\t').append(str(call))
                dosage = vcf_portion.vcf_matrix.get_matrix()[index_snp][index_donor]
                string_builder.append('\t').append(str(dosage))

            # add a newline
            string_builder.append('\n')
            # write the line to the file
            # TODO actually to file
            file.write(string_builder.to_string())

        # close the file we opened
        file.close()

        @staticmethod
        def get_dosages(vcf_portion, full_output_loc, ids_to_names=True):
            if vcf_portion.is_sorted is not True:
                vcf_portion.harmonize
            # we will write to a file
            string_builder = StringBuilder.StringBuilder()
            if ids_to_names:
                # get the header
                header = '\t'.join(vcf_portion.vcf_matrix.donor_ids)
                string_builder.append('\t').append(header).append('\n')
            else:
                # instead of the ids, write each fetched name instead
                for donor_id in vcf_portion.vcf_matrix.donor_ids:
                    donor_name = vcf_portion.vcf_matrix.donor_name_mapping[donor_id]
                    string_builder.append('\t').append(donor_name)

            # end the header
            string_builder.append('\n')
            # now each dose
            for snp_index in range(0, vcf_portion.vcf_matrix.snps.size):
                # we will write to a file, and build each line from a row in the matrix
                string_builder = StringBuilder.StringBuilder()
                # grab the name of the SNP
                if ids_to_names:
                    snp_name = vcf_portion.vcf_info.get_info('ID')[snp_index]
                    string_builder.append(snp_name)
                else:
                    snp_name = vcf_portion.vcf_info[snp_index]
                    string_builder.append(snp_name)
                # now get info for each participant
                line = '\t'.join(vcf_portion.vcf_matrix[snp_index])
                string_builder.append(line).append('\n')
                # write to the file







