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

    def set_category(self, row):
        """
        Adds runs to the tumors or hc lists, depending on whether the run is a tumor tissue or a healthy control
        :param row:         the row from the data frame
        :return:
        """
        # Make run objects
        one_run = Run(row['sample_num'], row['type_sample'], row['name'], row['category'])
        # Checked under which category the run falls and adds it to that specific list
        if row['category'] == 'tumor':
            self.tumors.append(one_run)
        elif row['category'] == 'hc':
            self.hc.append(one_run)


def prepair_file(path_file):
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


def make_objects(df_selection, dict_sample):
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


def main():
    # path_file = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/EGAD00001000292_metadata/delimited_maps/Sample_File.map"
    # head_path = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/"
    path_file = "D:/Hanze_Groningen/STAGE/non-coding-somatic-mutations-in-cancer/Anne/data/EGAD00001000292_Sample_File.map"
    df_selection = prepair_file(path_file)
    print(df_selection.head())
    dict_sample = dict()
    make_objects(df_selection, dict_sample)


if __name__ == "__main__":
    main()