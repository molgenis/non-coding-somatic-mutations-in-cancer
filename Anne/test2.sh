import pandas as pd
from itertools import combinations


class Run:
    """

    """

    def __init__(self, sample_num, type_sample, name, category):
        """
        Constructor
        :param sample_num:      sample number of the participant
        :param type_sample:     sample type (GL, FL, tFL)
        :param name:            specific number of the tissue
        :param category:        is the tissue tumor tissue or healthy tissue
        """
        self.sample_num = sample_num
        self.type_sample = type_sample
        self.name = name
        self.category = category


class Sample:
    """

    """

    def __init__(self):
        """
        Constructor
        """
        self.tumors = list()
        self.hc = list()
        self.sample_name = ''

    def set_category(self, row):
        """
        Adds runs to the tumors or hc lists, depending on whether the run is a tumor tissue or a healthy control
        :param row:         the row from the data frame
        :return:
        """
        # Make run objects
        one_run = Run(row['sample_num'], row['type_sample'], row['name'], row['category'])
        # Names sample_name
        self.sample_name = row['sample_num']
        # Checked under which category the run falls and adds it to that specific list
        if row['category'] == 'tumor':
            self.tumors.append(one_run)
        elif row['category'] == 'hc':
            self.hc.append(one_run)

    def get_arguments(self, head_path, type_sample='both'):
        """
        Causes all arguments to Mutect2/FilterMutectCalls to be made with all runs of a sample.
        :param head_path:   Main path to where the file will be saved
        :return:
        """
        # Starting string for file name/directory name
        arg = f'{self.sample_name}_'
        # Starting string for arg_mutect2 (the whole argument to eventually run Mutect2)
        arg_mutect2 = ''
        # Add each tumor to the different parameters
        for tum in self.tumors:
            #TODO split FL en tFL als gebruiker dat wil
            if type_sample='both':
                # Add tumor to arg (file name/directory name)
                arg += f'{tum.name}_'
                # For example, makes SS6005044 -> 5044
                number_tum = tum.name.replace("SS600", "")
                # Add tumor to arg_mutect2 (the whole argument to eventually run Mutect2)
                arg_mutect2 += f'-I {head_path}{number_tum}/SN_{tum.name}.DR.bam '
            elif type_sample='tFL':
                if tum.type_sample == 'tFL':
                    # Add tumor to arg (file name/directory name)
                    arg += f'{tum.name}_'
                    # For example, makes SS6005044 -> 5044
                    number_tum = tum.name.replace("SS600", "")
                    # Add tumor to arg_mutect2 (the whole argument to eventually run Mutect2)
                    arg_mutect2 += f'-I {head_path}{number_tum}/SN_{tum.name}.DR.bam '
            elif type_sample='FL':
                if tum.type_sample != 'tFL':
                    # Add tumor to arg (file name/directory name)
                    arg += f'{tum.name}_'
                    # For example, makes SS6005044 -> 5044
                    number_tum = tum.name.replace("SS600", "")
                    # Add tumor to arg_mutect2 (the whole argument to eventually run Mutect2)
                    arg_mutect2 += f'-I {head_path}{number_tum}/SN_{tum.name}.DR.bam '



        # Add each hc to the different parameters
        for hc in self.hc:
            # Add hc to arg (file name/directory name)
            arg += f'{hc.name}_'
            # For example, makes SS6005042 -> 5042
            number_hc = hc.name.replace("SS600", "")
            # Add hc to arg_mutect2 (the whole argument to eventually execute Mutect2)
            arg_mutect2 += f'-I {head_path}{number_hc}/SN_{hc.name}.DR.bam -normal {hc.name}.DR.bam '
        # Pastes everything together and makes three good arguments that can be passed to Mutect2/FilterMutectCalls
        arg_mutect2 += f'-O {head_path}{arg[:-1]}/unfiltered_{arg}somatic.vcf.gz'
        filter_mutect_calls_input = f'\n{head_path}{arg[:-1]}/unfiltered_{arg[:-1]}'
        filter_mutect_calls_output = f'\n{head_path}{arg[:-1]}/filtered_{arg[:-1]}'
        # File name (and path) after which these arguments are written
        name_file = f'{head_path}{arg[:-1]}.txt'
        # Calls the function that actually writes the arguments
        self.write_file(name_file, [arg_mutect2, filter_mutect_calls_input, filter_mutect_calls_output])

    def get_arguments_per_two(self, head_path, number_of_tumors=len(self.tumors), number_of_hc=len(self.hc)):
        """
        Ensures that all arguments for Mutect2/FilterMutectCalls are made with all possible combinations
        tumor vs hc of a sample (i.e. per two).
        :param head_path:   Main path to where the file will be saved
        :return:
        """
        number_of_tumors = check_length_list(self.tumors, number_of_tumors)
        number_of_hc = check_length_list(self.hc, number_of_hc)

        for edited_list_hc in combinations(self.hc, number_of_hc):
            for hc in edited_list_hc:
                # Add hc to arg (file name/directory name)
                arg += f'{hc.name}_'
                # For example, makes SS6005042 -> 5042
                number_hc = hc.name.replace("SS600", "")
                # Add hc to arg_mutect2 (the whole argument to eventually execute Mutect2)
                arg_mutect2 += f'-I {head_path}{number_hc}/SN_{hc.name}.DR.bam -normal {hc.name}.DR.bam '
            for edited_list_tumors in combinations(lst, number_of_tumors):
                for tum in edited_list_tumors:
                    #TODO split FL en tFL als gebruiker dat wil
                    if type_sample='both':
                        # Add tumor to arg (file name/directory name)
                        arg += f'{tum.name}_'
                        # For example, makes SS6005044 -> 5044
                        number_tum = tum.name.replace("SS600", "")
                        # Add tumor to arg_mutect2 (the whole argument to eventually run Mutect2)
                        arg_mutect2 += f'-I {head_path}{number_tum}/SN_{tum.name}.DR.bam '
                    elif type_sample='tFL':
                        if tum.type_sample == 'tFL':
                            # Add tumor to arg (file name/directory name)
                            arg += f'{tum.name}_'
                            # For example, makes SS6005044 -> 5044
                            number_tum = tum.name.replace("SS600", "")
                            # Add tumor to arg_mutect2 (the whole argument to eventually run Mutect2)
                            arg_mutect2 += f'-I {head_path}{number_tum}/SN_{tum.name}.DR.bam '
                    elif type_sample='FL':
                        if tum.type_sample != 'tFL':
                            # Add tumor to arg (file name/directory name)
                            arg += f'{tum.name}_'
                            # For example, makes SS6005044 -> 5044
                            number_tum = tum.name.replace("SS600", "")
                            # Add tumor to arg_mutect2 (the whole argument to eventually run Mutect2)
                            arg_mutect2 += f'-I {head_path}{number_tum}/SN_{tum.name}.DR.bam '

        
        for i in combinations(lst, lengthOfStrings):
            print(i)

    def check_length_list(self, number_list, number, type_sample):
        if type_sample == 'both':
            if len(number_list) < number:
                number = len(number_list)
        elif type_sample =='tFL':

        elif type_sample =='FL':
        
        return number_list, number

    def write_file(self, name_file, list_args):
        """
        Writes the arguments (list_args) to a file.
        :param name_file:   File name (and path) after which the arguments are written
        :param list_args:   list of arguments with the following arguments:
                            arg_mutect2:                    Argument to run Mutect2
                            filter_mutect_calls_input:      input file for FilterMutectCalls (after Mutect2)
                            filter_mutect_calls_output:     output file for FilterMutectCalls (after Mutect2)
        :return:
        """
        file_arguments = open(name_file, "w")
        # Loop over all arguments and writes these arguments line by line in a file
        for arg in list_args:
            file_arguments.write(f'{arg}')
        file_arguments.close()


def filter_file(path_file):
    """
    Creates a new dataframe with the information needed to create objects
    :param path_file:       the path to the file
    :return: df_selection:  the custom data frame from which to create objects and a dictionary
    """
    # Read file
    df = pd.read_csv(path_file, sep='\t', header=None)
    # Make empty dataframe
    df_selection = pd.DataFrame()
    # Split first column on _ in sample_num and type_sample (for example S1_FL --> S1   FL)
    df_selection[['sample_num', 'type_sample']] = df[0].str.split("_", expand=True, )
    # Take only the number of the file (eg SS6005044)
    df_selection['name'] = df[2].str.split(".", expand=True, )[0]
    # Make extra column with whether something is a tumor or hc (healthy control)
    df_selection.loc[df_selection['type_sample'].str.contains('FL'), 'category'] = 'tumor'
    df_selection.loc[df_selection['type_sample'].str.contains('GL'), 'category'] = 'hc'
    return df_selection


def make_objects(df_selection, head_path):
    """
    Creates the dictionary with the objects per participant/sample
    :param df_selection:    the custom data frame from which to create objects and a dictionary
    :param head_path:       Main path to where the file will be saved
    :return:                a dictionary containing the sample_num as key and as value Sample objects
    """
    # Empty dict:
    dict_samples = dict()
    # Loop over lines in data frame
    for index, row in df_selection.iterrows():
        # Checked if a key is already in the dictionary or not
        if row['sample_num'] in dict_samples.keys():
            # update object Sample
            dict_samples[row['sample_num']].set_category(row)
        else:
            # Make object Sample
            one_sample = Sample()
            one_sample.set_category(row)
            dict_samples[row['sample_num']] = one_sample
    return dict_samples


def arguments_to_file(dict_samples, head_path):
    """
    Ensures that all arguments are created and written.
    :param dict_samples:     a dictionary containing the sample_num as key and as value Sample objects
    :param head_path:       Main path to where the file will be saved
    :return:
    """
    # Loop over keys from dict_sample
    for key in dict_samples:
        dict_samples[key].get_arguments(head_path)
        dict_samples[key].get_arguments_per_two(head_path)


def main():
    """

    :return:
    """
    # the path to the file
    # path_file = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/EGAD00001000292_metadata/delimited_maps/Sample_File.map"
    # Main path to where the file will be saved
    # head_path = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/"
    # the path to the file
    path_file = "D:/Hanze_Groningen/STAGE/non-coding-somatic-mutations-in-cancer/Anne/data/EGAD00001000292_Sample_File.map"
    # Main path to where the file will be saved
    head_path = "D:/Hanze_Groningen/STAGE/non-coding-somatic-mutations-in-cancer/Anne/data/"
    # Filter file
    df_selection = filter_file(path_file)
    print(df_selection)
    # Make objects for all samples/participants
    #dict_samples = make_objects(df_selection, head_path)
    #arguments_to_file(dict_samples, head_path)


if __name__ == "__main__":
    main()
