'''
Created on 12 feb. 2022

@author: roy.oelen
'''

# external imports
import multiprocessing
import random
import string
# internal imports
from multiprocessing import Pool

import DbConnector
import FileWriter
import TypeConverter


class JobCluster:
    '''
    object that is an abstraction layer for the communication with the database and handling the returns
    '''


    def __init__(self, database_loc, num_threads=1, num_db_copies=1, is_test=False):
        '''
        constructor
        :param database_loc: base path of the database (before '.db')
        :param n_copies: the number of copies of the databaase are available
        :param is_test: convenience key to make testing easier
        '''
        self.num_threads = num_threads
        self.num_db_copies = num_db_copies
        self.base_connection = None
        self.is_test = is_test
        self.database_loc = database_loc
        self.setup_cluster(database_loc, num_db_copies)

    def setup_cluster(self, database_loc, n_copies):
        '''
        set up the number of database file copies and paths
        :param database_loc: base path of the database (before '.db')
        :param n_copies: the number of copies of the databaase are available
        :return: NONE
        '''
        # check if anything set
        if n_copies is None:
            raise Exception('number of connections cannot be None')
        elif isinstance(n_copies, int) is False:
            raise Exception('number of connections must be numeric')
        elif n_copies < 0:
            raise Exception('number of connections cannot a negative number')
        elif n_copies == 0:
            raise Exception('there must be at least one connection')
        elif n_copies >= 1:
            self.num_threads = n_copies
            # set up the base connection for lightweight queries
            database_loc_lightweight = ''.join([self.database_loc, '.db'])
            base_connection = DbConnector.DbConnector(database_loc_lightweight)
            self.base_connection = base_connection

    def get_studies(self):
        '''
        get the studies available
        :return: a tuple of the internal and external IDs of each study
        '''
        # we will just use the first connection
        studies = self.base_connection.get_studies()
        return studies

    def show_tables(self):
        '''
        get the tables in the database
        :return: a tuple with the table and their properties
        '''
        # just the first connection
        tables = self.base_connection.get_tables()
        return tables

    def get_columns(self, table_name):
        '''
        get the columns for a specific table
        :param table_name: the table you want the data of
        :return: a tuple describing the columns and properties of the table
        '''
        columns = self.base_connection.get_columns(table_name)
        return columns

    def database_to_vcf(self, study_id, full_output_loc, alternate_database_loc=None):
        '''
        convert the data in the database of a specific study to a VCF file
        :param study_id: the internal ID of the study you want to create a VCF for
        :param full_output_loc: the base output loc of the resulting VCF
        :return:
        '''
        # connect to the database
        print(''.join(['starting study', str(study_id)]))
        db_connector = None
        if alternate_database_loc is not None:
            db_connector = DbConnector.DbConnector(alternate_database_loc)
        else:
            db_connector = DbConnector.DbConnector(self.database_loc)
        # get the results for the study
        study_result = db_connector.get_study_data(study_id)
        print(study_result)
        FileWriter.FileWriter.write_vcf(study_result, full_output_loc)
        return True

    @staticmethod
    def get_study_data_static(study_id, full_output_loc, database_loc):
        # connect to the database
        print(''.join(['starting study', str(study_id)]))
        db_connector = DbConnector.DbConnector(database_loc)
        # get the results for the study
        study_result = db_connector.get_study_data(study_id)
        return study_result

    @staticmethod
    def database_to_vcf_static(study_id, full_output_loc, database_loc):
        # get the study data
        study_result = JobCluster.get_study_data_static(study_id, full_output_loc, database_loc)
        # harmonize the object so that info and db are in the same order
        study_result.harmonize()
        # start writing to vcf files
        TypeConverter.TypeConverter.vcf_portion_to_file(study_result, full_output_loc)



    def studies_to_vcf(self, base_output_loc, studies=[]):
        '''
        convert the data in the database of a set of studies to a VCF file
        :param base_output_loc: the base output loc of the resulting VCF
        :param studies: a lisft of the internal IDs of the studies you want to write to VCF
        :return: NONE
        '''
        # try the default first
        studies_to_use = studies
        # if there are no studies, we'll just do all
        if len(studies_to_use) < 1:
            studies_to_use = self.get_studies()
        # we will use a number of copies of the database
        copy_number = 1
        # let's get those jobs then
        jobs = []
        file_names = []
        for study in studies_to_use:
            # build the output location
            random_string = JobCluster.get_random_string(15)
            # paste the output together
            output_file_location = ''.join([base_output_loc, '/', random_string, '.vcf.tmp'])
            file_names.append(output_file_location)

            # get a database file to use
            database_copy_location = ''.join([self.database_loc, '_', str(copy_number), '.db'])
            # increment the counter
            copy_number = copy_number + 1
            # reset the count if required
            if copy_number > self.num_threads:
                copy_number = 1

            # make a process
            process = multiprocessing.Process(target=self.database_to_vcf, args=(study['study_id'], output_file_location,database_copy_location))
            # add to list of processes
            jobs.append(process)
            # finally start the job
            process.start()

        # now wait for jobs to be done
        for job in jobs:
            job.join()
        # merge everything
        self.merge_files(file_names)

    def studies_to_vcf_pooled(self, base_output_loc, studies=[], threads=multiprocessing.cpu_count()-1):
        '''
        convert the data in the database of a set of studies to a VCF file
        :param base_output_loc: the base output loc of the resulting VCF
        :param studies: a lisft of the internal IDs of the studies you want to write to VCF
        :return: NONE
        '''
        # try the default first
        studies_to_use = studies
        # if there are no studies, we'll just do all
        if len(studies_to_use) < 1:
            studies_to_use = self.get_studies()

        print(studies_to_use)
        # we will use a number of copies of the database
        copy_number = 1
        # make a pool
        pool = Pool(processes=threads)
        # check each study
        for study in studies_to_use:
            # build the output location
            random_string = JobCluster.get_random_string(15)
            # paste the output together
            output_file_location = ''.join([base_output_loc, '/', random_string, '.vcf.tmp'])

            # get a database file to use
            database_copy_location = ''.join([self.database_loc, '_', str(copy_number), '.db'])
            # increment the counter
            copy_number = copy_number + 1
            # reset the count if required
            if copy_number > self.num_db_copies:
                copy_number = 1

            # add to the pool
            #pool.apply_async(JobCluster.database_to_vcf_static, args=(study['study_id'], output_file_location, database_copy_location,))
            pool.apply(JobCluster.database_to_vcf_static, args=(study['study_id'], output_file_location, database_copy_location))


        pool.close()
        pool.join()

    def merge_files(self, file_names):
        '''

        :param file_names:
        :return:
        '''
        # TODO actually do this
        print('merge_files not yet implemented')
        None

    @staticmethod
    def get_random_string(string_length):
        '''
        get a random string of ascii letters
        :param string_length: how long the random string should be
        :return: a random string of ascii letters of the supplied length
        '''
        # use random to get a random string
        random_string = ''.join(random.choice(string.ascii_letters) for i in range(string_length))
        return random_string
