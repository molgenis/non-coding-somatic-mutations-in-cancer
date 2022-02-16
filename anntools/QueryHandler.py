'''
Created on 12 feb. 2022

@author: roy.oelen
'''

# external imports
import argparse
# internal imports
import datetime

import JobCluster


class QueryHandler:
    '''
    handle requests from the database to other formats
    '''

    def __init__(self, db_location, output_location, study, num_threads, cache_loc, is_test=False):
        '''
        contructor
        :param db_location: location of the database file
        :param output_location: where to write the resulting output file
        :param study: the study you want to grab data from
        :param num_threads: the max number of threads to use
        :param cache_loc: location where the temporary files are stored before merging od data
        '''
        # initialize the object with relevant
        self.db_location = db_location
        self.output_location = output_location
        self.study = study
        self.num_threads = num_threads
        self.cache_loc = cache_loc
        self.job_manager = None
        self.is_test = is_test
        self.study_cache = {}
        self.setup_job_manager()

    def setup_job_manager(self):
        '''
        set up the cluster of database connections

        :return: NONE does not return anything
        '''
        self.job_manager = JobCluster.JobCluster(self.db_location, self.num_threads, is_test=self.is_test)

    def get_studies(self):
        '''
        get the available studies
        :return: tuples of the studies, with the internal and external identifier
        '''
        # we can keep track of the studies without going through the database each time
        if(len(self.study_cache.keys()) == 0):
            # however the first time we need to set up the cache
            self.study_cache = self.job_manager.get_studies()

        # regardless of whether it was filled or not, it should be filled now
        studies = self.study_cache
        return studies

    def clear_study_cache(self):
        '''
        empty the cache of the study tuples (if changes have been made since starting this application)
        :return:
        '''
        # for whatever reason, we need to clear cache
        self.study_cache = {}

    def show_tables(self):
        '''
        show the tables present in the database
        :return: the tables present in the database
        '''
        tables = self.job_manager.show_tables()
        return tables


    def get_columns(self, table_name):
        '''
        get the columns present for a specific table
        :param table_name: the table to check the columns for
        :return:
        '''
        column_names = self.job_manager.get_columns(table_name)
        return column_names

    def studies_to_vcf(self, base_output_loc, studies=[]):
        '''
        convert studies in the database to vcf files
        :param base_output_loc: the base location of the database files
        :param studies: the studies to generate VCF files for
        :return: NONE a resulting vcf file is made at the required location
        '''
        self.job_manager.studies_to_vcf(base_output_loc, studies)






def construct_query_handler():
    '''
    create a QueryHandler to abstract the database actions

    :return: a QueryHandler object
    '''
    args = parse_args()
    query_handler = QueryHandler(
        db_location=args.db_location,
        output_location=args.output_location,
        study = args.study,
        num_threads = args.num_threads,
        cache = args.cache_loc
    )


def parse_args():
    '''
    parse the command line arguments

    :return: a dictionary with keys and values like they would be created by argparse
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--db-location', type = str, help = 'location of the sqlite database (string)')
    parser.add_argument('-o', '--output-location', type = str, help = 'location to store the output (string)')
    parser.add_argument('-s', '--study', type = str, help = 'name of the study to fetch (string), default=*', default = '*')
    parser.add_argument('-t', '--num-threads', type = int, help = 'number of threads to use (numeric), default=1', default=1)
    parser.add_argument('-c', '--cache-loc', type = str, help = 'location for caching temporary files (string), default=~/', default='~/')
    parser.parse_args()
    return parser


def test_run():
    '''
    test the application with standard paths
    :return: NONE
    '''

    # create a dictionary that looks like what the argparse would look like
    setup = {'db_location' : '/Users/royoelen/Desktop/Database_internship_UPDATE',
             'output_location' : '/Users/royoelen/Desktop/',
             'study' : 'bla',
             'num_threads' : 2,
             'cache_loc' : '/Users/royoelen/Desktop/'}
    # create the object
    query_handler = QueryHandler(
        db_location=setup['db_location'],
        output_location=setup['output_location'],
        study = setup['study'],
        num_threads = setup['num_threads'],
        cache_loc = setup['cache_loc'],
        is_test = True
    )
    # use a couple of functions
    #query_handler.show_tables()
    #query_handler.get_studies()
    #query_handler.get_columns('project')
    #query_handler.studies_to_vcf('/Users/royoelen/Desktop/test.vcf',
    #                               [{'study_id' : 1, 'study_name' : 'ALL-US'},
    #                               {'study_id' : 2, 'study_name' : 'AML-US'}])

    query_handler.studies_to_vcf('/Users/royoelen/Desktop/test.vcf')

if __name__ == '__main__':
    '''
    application start
    '''
    test_run()
    #construct_query_handler()




