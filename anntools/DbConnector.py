'''
Created on 12 feb. 2022

@author: roy.oelen
'''

# own imports
import TypeConverter
import Constants
import VcfInfo
import VcfMatrix
import VcfPortion

# external exports
import sqlite3


class DbConnector:
    '''
    Connect to the database and fetch data
    '''

    def __init__(self, database_loc=None):
        '''
        Constructor
        '''
        self.connection = None
        self.cursor = None
        self.database_loc = database_loc
        # connect database
        if self.database_loc is not None:
            self.connect()

    def connect(self):
        '''
        connect to a database
        '''
        self.connection = sqlite3.connect(self.database_loc)
        self.cursor = self.connection.cursor()

    def get_studies(self):
        '''
        get the the studies in the database, with the internal and external IDs
        :return:
        '''
        # get the correct table name
        study_table_name = Constants.NAME_TABLE_STUDIES
        studies_table = DbConnector.sanitize_string(study_table_name)
        # create a query
        query = ''.join(['SELECT * FROM ', studies_table])
        # we will have an entry for each study
        studies=[]
        # execute the query
        self.cursor.execute(query)
        # change each row into an entry in our list
        for row in self.cursor.fetchall():
            study = {'study_id' : row[0], 'study_name' : row[1]}
            studies.append(study)
        return studies


    def get_tables(self):
        '''
        get the tables in the database
        :return: a tuple with the table and their properties
        '''
        query = 'SELECT name FROM sqlite_schema WHERE type = \'table\' AND name NOT LIKE \'sqlite+%\''
        # get all the tables
        for table in self.cursor.execute(query):
            print(table)

    def get_columns(self, table_name):
        '''
        get the columns for a specific table
        :param table_name: the table you want the data of
        :return: a tuple describing the columns and properties of the table
        '''
        sanitized_table_name = DbConnector.sanitize_string(table_name)
        query = 'PRAGMA table_info({})'
        self.cursor.execute(query.format(sanitized_table_name))
        rows = self.cursor.fetchall()
        return rows

    @staticmethod
    def sanitize_string(string):
        '''
        remove illegal characters
        :param string:
        :return:
        '''
        # TODO actually sanitize...
        sanitized = string
        return sanitized

    def get_study_data(self, study_id):
        '''
        get the data of one study into a vcf compatible object
        :param study_id: the internal identifier of the study you want the data of
        :return:
        '''
        # get the donors for this study
        donors = self.get_donors_from_study(study_id)
        # get just the raw participant id
        donor_ids = []
        for entry in donors:
            donor_ids.append(entry[0])
        # get the SNP metadata
        snp_id_tuples = self.get_snp_ids_for_donor_ids(donor_ids)
        snp_ids = []
        # the raw ids again
        for entry in snp_id_tuples:
            snp_ids.append(entry[0])
        # fetch the metadata for these snps
        snp_metadatas = self.get_snp_metadata(snp_ids)
        # convert to a VcfInfo object
        vcf_info = TypeConverter.TypeConverter.tuples_to_vcfinfo(snp_metadatas)
        # now get the genotype data
        snp_genodata = self.get_snp_genodata(snp_ids, donor_ids)
        # turn into VcfMatrix object
        vcf_geno = TypeConverter.TypeConverter.tuples_to_vcfdata(snp_genodata, snp_ids, donor_ids)
        # turn into a VcfPortion object
        vcf_portion = VcfPortion.VcfPortion(vcf_geno, vcf_info)
        # for memory sake, let's remove some variables
        del(snp_metadatas)
        del(snp_genodata)
        # return the result
        return vcf_portion




    def get_donors_from_study(self, study_id):
        '''
        get the internal and external ids of all the donors included in a study
        :param study_id: the internal id of the study you want to get the donors for
        :return: the internal and external ids of all the donors included
        '''
        # create a query to get the donors
        query = 'SELECT d.ID, d.donor_ID FROM donor d WHERE d.project_ID=?'
        #query = 'SELECT \'d.ID\', \'d.donor_ID\', FROM \'donor\' d WHERE \'d.project_ID\'=?'
        # add the project ID we want
        self.cursor.execute(query, (int(study_id),))
        # fetch the result
        rows = self.cursor.fetchall()
        # return the result
        return rows

    def get_snp_ids_for_donor_ids(self, donor_ids):
        '''
        get IDs of all the snps for a list of donors
        :param donor_ids: a list with the IDs of the donors we want to get the (outer join of the) SNP ids of
        :return: a tuple of the snp IDs
        '''
        # create a query
        query_start = 'SELECT DISTINCT dh.snp_ID FROM donor_has_snp dh WHERE dh.donor_ID IN ('
        # add the donor ID questions
        query_in_parameters = ','.join('?' * len(donor_ids))
        # finish the query
        query = ''.join([query_start, query_in_parameters, ')'])
        # execute the query
        self.cursor.execute(query, tuple(donor_ids))
        # fetch the result
        rows = self.cursor.fetchall()
        # return the result
        return rows

    def get_snp_metadata_donors(self, donor_ids):
        '''
        get the meta data for the snps of the donors
        :param donor_ids: a list of donors we want the (outer join of the) SNPs metadata of
        :return: a typle of the metadata for each SNP
        '''
        # create a query to select
        query_start = 'SELECT snp.ID, snp.chr, snp.pos_start, snp.pos_end, snp.ref, snp.alt FROM snp snp JOIN donor_has_snp dhs ON snp.ID = dhs.snp_ID WHERE dhs.snp_ID IN ('
        # add the donor ID questions
        query_in_parameters = ','.join('?' * len(donor_ids))
        # finish the query
        query = ''.join([query_start, query_in_parameters, ')'])
        # execute the query
        self.cursor.execute(query, tuple(donor_ids))
        # fetch the result
        rows = self.cursor.fetchall()
        return rows

    def get_snp_metadata(self, snp_ids):
        '''
        get the meta data for the snps supplied
        :param snp_ids: a list of the internal IDs of the SNPs
        :return: a typle of the metadata for each SNP
        '''
        # create a query to select
        query_start = 'SELECT snp.ID, snp.chr, snp.pos_start, snp.pos_end, snp.ref, snp.alt FROM snp snp WHERE snp.ID IN ('
        # add the donor ID questions
        query_in_parameters = ','.join('?' * len(snp_ids))
        # finish the query
        query = ''.join([query_start, query_in_parameters, ')'])
        # execute the query
        rows = None
        # check the size of the query
        if len(snp_ids) > 999:
            # we need to do this multple times then
            list_results = []
            # we will start at zero
            start = 0
            end = 0
            while end < len(snp_ids):
                # check how far to go
                end = start + 999
                # of course we won't end in a exactly the right sizes of slices
                if end > len(snp_ids):
                    # we will grab to the end of the list, instead of start + 1000
                    end = start + len(snp_ids)
                # get the relevant ids
                ids_interation = snp_ids[start:end]
                # use these for our query
                query_in_parameters = ','.join('?' * len(ids_interation))
                # finish the query
                query = ''.join([query_start, query_in_parameters, ')'])
                self.cursor.execute(query, tuple(ids_interation))
                # fetch the result
                rows_portion = self.cursor.fetchall()
                # add to the list of queries
                list_results.append(rows_portion)
                # the new starts is the old start
                start = end + 1

            # merge all of the lists of lists, into one list
            rows = [j for i in list_results for j in i]

        else:
            self.cursor.execute(query, tuple(snp_ids))
            # fetch the result
            rows = self.cursor.fetchall()
        return rows

    def get_snp_genodata(self, snp_ids, donor_ids, max_snp_id_length=499, max_donor_id_length=499):
        # init the start of the query
        query_start = 'SELECT dhs.snp_ID, dhs.donor_ID, dhs.mutant_allele_read_count, dhs.total_read_count FROM donor_has_snp dhs WHERE dhs.snp_ID IN ('
        # add the donor ID questions
        query_in_parameters = ','.join('?' * len(snp_ids))
        # finish the query
        query = ''.join([query_start, query_in_parameters, ') AND donor_ID IN ('])
        # execute the query
        rows = None
        # check the size of the query
        if len(snp_ids) > 999:
            # we need to do this multple times then
            list_results = []
            # we will start at zero
            start = 0
            end = 0
            while end < len(snp_ids):
                # init because we can split a second time
                rows_portion = None
                # check how far to go
                end = start + max_snp_id_length
                # of course we won't end in a exactly the right sizes of slices
                if end > len(snp_ids):
                    # we will grab to the end of the list, instead of start + 1000
                    end = start + len(snp_ids)
                # get the relevant ids
                ids_interation = snp_ids[start:end]
                # use these for our query
                query_in_parameters = ','.join('?' * len(ids_interation))
                # finish the query
                query_first_part = ''.join([query_start, query_in_parameters, ') AND donor_ID IN ('])
                # check if we also need to subset the donors
                if len(donor_ids) > max_donor_id_length:
                    rows_portion = self.split_query(query_first_part, ids_interation, donor_ids, max_donor_id_length)
                else:
                    # build the rest of the query
                    query_in_parameters_2 = ','.join('?' * len(donor_ids))
                    # finish the query
                    query = ''.join([query_first_part, query_in_parameters_2, ')'])
                    self.cursor.execute(query, tuple(ids_interation + donor_ids))
                    # fetch the result
                    rows_portion = self.cursor.fetchall()

                # add to the list of queries
                list_results.append(rows_portion)
                # the new starts is the old start
                start = end + 1
            # merge all of the lists of lists, into one list
            rows = [j for i in list_results for j in i]
        else:
            # we don't need to split by snps, but we could still need to split the donors
            if len(donor_ids) > max_donor_id_length:
                rows = self.split_query(query, donor_ids, snp_ids, max_donor_id_length)
            else:
                # build the rest of the query
                query_in_parameters_2 = ','.join('?' * len(donor_ids))
                # finish the query
                query = ''.join([query, query_in_parameters_2, ')'])
                self.cursor.execute(query, tuple(snp_ids + donor_ids))
                # fetch the result
                rows = self.cursor.fetchall()

        return rows

    def split_query(self, query_start, ids, query_start_ids, max_id_length=499):
        # we need to do this multiple times then
        list_results = []
        # we will start at zero
        start = 0
        end = 0
        while end < len(ids):
            # check how far to go
            end = start + max_id_length
            # of course we won't end in a exactly the right sizes of slices
            if end > len(ids):
                # we will grab to the end of the list, instead of start + 1000
                end = start + len(ids)
            # get the relevant ids
            ids_interation = ids[start:end]
            # use these for our query
            query_in_parameters = ','.join('?' * len(ids_interation))
            # finish the query
            query = ''.join([query_start, query_in_parameters, ')'])
            self.cursor.execute(query, tuple(query_start_ids, ids_interation))
            # fetch the result
            rows_portion = self.cursor.fetchall()
            # add to the list of queries
            list_results.append(rows_portion)
            # the new starts is the old start
            start = end + 1
        rows = [j for i in list_results for j in i]
        return rows


    def create_temporary_donor_table(self, study_id):
        '''

        :param study_id:
        :return:
        '''
        # TODO finished method
        print('create_temporary_donor_table might fail')

        # make sure the input is clean
        study_id_to_use = self.sanitize_string(str(study_id))
        # create a temporary table first
        query = 'CREATE TEMP TABLE donor_' + study_id_to_use + ' AS ' \
                                                               'SELECT * ' \
                                                               'FROM donor ' \
                                                               'WHERE project_ID=' + study_id_to_use+ ';'
        print(query)
        #self.cursor.execute(query, (study_id_to_use))
        self.cursor.execute(query)
        self.connection.commit()
        print(self.get_tables())

