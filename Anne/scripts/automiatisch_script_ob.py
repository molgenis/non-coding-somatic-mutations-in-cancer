import pandas as pd


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
        self.sample_name = row['sample_num']
        # Checked under which category the run falls and adds it to that specific list
        if row['category'] == 'tumor':
            self.tumors.append(one_run)
        elif row['category'] == 'hc':
            self.hc.append(one_run)

    def get_arguments(self, head_path):
        """

        :param head_path:
        :return:
        """
        arg = f'{self.sample_name}_'
        arg_mutect2 = ''
        for tum in self.tumors:
            arg += f'{tum.name}_'
            number_tum = tum.name.replace("SS600", "")
            arg_mutect2 += f'-I {head_path}{number_tum}/SN_{tum.name}.DR.bam '
        for hc in self.hc:
            arg += f'{hc.name}_'
            number_hc = hc.name.replace("SS600", "")
            arg_mutect2 += f'-I {head_path}{number_hc}/SN_{hc.name}.DR.bam -normal {hc.name}.DR.bam '

        arg_mutect2 += f'-O {head_path}{arg[:-1]}/unfiltered_{arg}somatic.vcf.gz'
        filter_mutect_calls_input = f'\n{head_path}{arg[:-1]}/unfiltered_{arg[:-1]}'
        filter_mutect_calls_output = f'\n{head_path}{arg[:-1]}/filtered_{arg[:-1]}'
        name_file = f'{head_path}{arg[:-1]}.txt'
        self.write_file(name_file, [arg_mutect2, filter_mutect_calls_input, filter_mutect_calls_output])

    def get_arguments_per_two(self, head_path):
        """

        :param head_path:
        :return:
        """
        for hc in self.hc:
            number_hc = hc.name.replace("SS600", "")
            for tum in self.tumors:
                number_tum = tum.name.replace("SS600", "")
                arg = f'{self.sample_name}_{tum.name}_{hc.name}'
                arg_mutect2 = f'-I {head_path}{number_tum}/SN_{tum.name}.DR.bam -I {head_path}{number_hc}/SN_{hc.name}.DR.bam -normal {hc.name}.DR.bam -O {head_path}{arg}/unfiltered_{arg}_somatic.vcf.gz'
                filter_mutect_calls_input = f'\n{head_path}{arg}/unfiltered_{arg}'
                filter_mutect_calls_output = f'\n{head_path}{arg}/filtered_{arg}'
                name_file = f'{head_path}{arg}.txt'
                self.write_file(name_file, [arg_mutect2, filter_mutect_calls_input, filter_mutect_calls_output])

    def write_file(self, name_file, list_args):
        """

        :param name_file:
        :param list_args:
        :return:
        """
        f = open(name_file, "w")
        for arg in list_args:
            f.write(f'{arg}')
        f.close()


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


def make_objects(df_selection, dict_sample, head_path):
    """
    Creates the dictionary with the objects per participant/sample
    :param df_selection:    the custom data frame from which to create objects and a dictionary
    :param dict_sample:     a dictionary containing the sample_num as key and as value Sample objects
    :return:
    """
    # Loop over lines in data frame
    for index, row in df_selection.iterrows():
        # Checked if a key is already in the dictionary or not
        if row['sample_num'] in dict_sample.keys():
            # update object Sample
            dict_sample[row['sample_num']].set_category(row)
        else:
            # Make object Sample
            one_sample = Sample()
            one_sample.set_category(row)
            dict_sample[row['sample_num']] = one_sample

    # for i in dict_sample['S2'].tumors:
    #     print(f'{i.sample_num} - {i.name}')
    # for i in dict_sample['S2'].hc:
    #     print(f'{i.sample_num} - {i.name}')

    #for key, value in d.items():
    for key in dict_sample:
        dict_sample[key].get_arguments(head_path)
        dict_sample[key].get_arguments_per_two(head_path)


def main():
    """
    
    :return:
    """
    #path_file = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/EGAD00001000292_metadata/delimited_maps/Sample_File.map"
    #head_path = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/"
    path_file = "D:/Hanze_Groningen/STAGE/non-coding-somatic-mutations-in-cancer/Anne/data/EGAD00001000292_Sample_File.map"
    head_path = "D:/Hanze_Groningen/STAGE/non-coding-somatic-mutations-in-cancer/Anne/data/"
    df_selection = filter_file(path_file)
    print(df_selection.head())
    dict_sample = dict()
    make_objects(df_selection, dict_sample, head_path)


if __name__ == "__main__":
    main()
