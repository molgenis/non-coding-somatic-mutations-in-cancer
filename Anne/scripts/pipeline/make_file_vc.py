import pandas as pd
from itertools import combinations
import sys


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

    def get_arguments(self, head_path, chrom, compare_hc_tum, manual_comparison, mutect2_comparison, number_of_tumors=None,
                      number_of_hc=None,
                      type_sample='both', type_aln='bowtie', type_aln2='bowtie2'):
        """
        Ensures that all arguments for Mutect2/FilterMutectCalls are made with all possible combinations
        tumor vs hc of a sample (i.e. per two).
        :param head_path:           Main path to where the file will be saved
        :param compare_hc_tum:
        :param manual_comparison:
        :param mutect2_comparison:
        :param number_of_tumors:    Number of hc you want to combine while running Mutect2
        :param number_of_hc:        Number of hc you want to combine while running Mutect2
        :param type_sample:         What type of tumor you want to have ("both", "tFL" or "FL")
        :param type_aln:
        :param type_aln2:
        :return:
        """
        if number_of_tumors is None:
            number_of_tumors = len(self.tumors)
        if number_of_hc is None:
            number_of_hc = len(self.hc)
        # Check if the number is possible.
        edited_list_tumors, number_of_tumors = self.check_length_list(self.tumors, number_of_tumors, type_sample)
        edited_list_hc, number_of_hc = self.check_length_list(self.hc, number_of_hc, type_sample, 'hc')

        for edited_number_hc in combinations(edited_list_hc, number_of_hc):
            # Starting string for file name/directory name
            arg_hc = f'{self.sample_name}_numT_{number_of_tumors}_numHC_{number_of_hc}_'
            # Name the chosen tumor type in the args. (so you can see the difference in files)
            if type_sample != 'both':
                arg_hc += f'{type_sample}_'
            # Starting string for arg_mutect2 (the whole argument to eventually run Mutect2)
            arg_mutect2_hc = ''
            # Add each hc to the different parameters
            for hc in edited_number_hc:
                # Add hc to arg (file name/directory name)
                arg_hc += f'{hc.name}_'
                # For example, makes SS6005042 -> 5042
                number_hc = hc.name.replace("SS600", "")
                # single mutect2
                self.write_file(f'{head_path}{self.sample_name}/{chrom}/mutect_{type_aln}/{type_aln}_{hc.name}.txt',
                                [f'{head_path}{self.sample_name}/{number_hc}_vcf/{chrom}/{type_aln}\n',
                                 f'-I {head_path}{self.sample_name}/{number_hc}/{chrom}/{type_aln}/SN_{hc.name}.bam '
                                 f'-normal {type_aln2}_{hc.name}.DR ',
                                 f'\n{head_path}{self.sample_name}/{number_hc}_vcf/{chrom}/{type_aln}/{hc.name}_'])
                # Add hc to arg_mutect2 (the whole argument to eventually execute Mutect2)
                arg_mutect2_hc += f'-I {head_path}{self.sample_name}/{number_hc}/{chrom}/{type_aln}/SN_{hc.name}.bam ' \
                                  f'-normal {type_aln2}_{hc.name}.DR '
            for edited_number_tumors in combinations(edited_list_tumors, number_of_tumors):
                arg_tumor = ''
                arg_mutect2_tumor = ''
                for tum in edited_number_tumors:
                    # Add tumor to arg (file name/directory name)
                    arg_tumor += f'{tum.name}_'
                    # For example, makes SS6005044 -> 5044
                    number_tum = tum.name.replace("SS600", "")
                    if number_of_hc == 1 and number_of_tumors == 1:
                        compare_hc_tum.write(
                            f"{head_path}{self.sample_name}/{number_hc}_vcf/{chrom}/{type_aln}/"
                            f"{hc.name}_somatic_filtered.vcf.gz "
                            f"{head_path}{self.sample_name}/{number_tum}_vcf/{chrom}/{type_aln}/"
                            f"{tum.name}_somatic_filtered.vcf.gz"
                            f" -p {head_path}{self.sample_name}/compare_{number_hc}_{number_tum}/{chrom}/{type_aln}/\n")
                        manual_comparison.write(
                            f"{head_path}{self.sample_name}/compare_{number_hc}_{number_tum}/{chrom}/{type_aln}/0001.vcf.gz "
                        )
                    # single mutect2
                    self.write_file(f'{head_path}{self.sample_name}/{chrom}/mutect_{type_aln}/{type_aln}_{tum.name}.txt',
                                    [f'{head_path}{self.sample_name}/{number_tum}_vcf/{chrom}/{type_aln}\n',
                                     f'-I {head_path}{self.sample_name}/{number_tum}/{chrom}/{type_aln}/SN_{tum.name}.bam ',
                                     f'\n{head_path}{self.sample_name}/{number_tum}_vcf/{chrom}/{type_aln}/{tum.name}_'])
                    # Add tumor to arg_mutect2 (the whole argument to eventually run Mutect2)
                    arg_mutect2_tumor += f'-I {head_path}{self.sample_name}/{number_tum}/{chrom}/{type_aln}/SN_{tum.name}.bam '
                # Merge hc and tumor arguments
                arg = arg_hc + arg_tumor
                arg_mutect2 = arg_mutect2_hc + arg_mutect2_tumor
                # Pastes everything together and makes three good arguments that can be passed
                # to Mutect2/FilterMutectCalls
                # Path for mkdir
                mkdir_path = f'{head_path}{self.sample_name}/{arg[:-1]}/{chrom}/{type_aln}\n'
                mutect2_comparison.write(f'{head_path}{self.sample_name}/{arg[:-1]}/{chrom}/{type_aln}/'
                                         f'{arg}somatic_filtered_PON_GERM.vcf.gz ')
                # Path and first part of output files
                file_output = f'\n{head_path}{self.sample_name}/{arg[:-1]}/{chrom}/{type_aln}/{arg}'
                # File name (and path) after which these arguments are written
                name_file = f'{head_path}{self.sample_name}/{chrom}/mutect_{type_aln}/{type_aln}_{arg[:-1]}.txt'

                # Calls the function that actually writes the arguments
                self.write_file(name_file, [mkdir_path, arg_mutect2, file_output])

    def check_length_list(self, number_list, number, type_sample, category=''):
        """
        Check if the number is possible. If not (i.e. if it is too large a number) then take the length of the list.
        :param number_list: List of tumors/hc
        :param number:      Number of tumors/hc you want to combine while running Mutect2
        :param type_sample: What type of tumor you want to have ("both", "tFL" or "FL")
        :param category:    Whether it is hc or tumor
        :return:            edit_list or number_list, whose edit_list is modified if a specific tumor type is chosen.
                            number: possibly adjusted number to the length of the list. (Number of tumors/hc you want to
                            combine while running Mutect2)
        """
        # Check which tumor type is chosen.
        # If it is both and if it is not tumor type but hc the following is performed.
        if type_sample == 'both' or category == 'hc':
            # Check if the number is possible
            if len(number_list) < number:
                number = len(number_list)
            return number_list, number
        # When a tumor type has been selected, the tumors are checked for type.
        # They are placed in a list of selected tumor types.
        else:
            edit_list = list()
            for value in number_list:
                if type_sample == 'tFL' and value.type_sample == 'tFL':
                    edit_list.append(value)
                elif type_sample == 'FL' and value.type_sample != 'tFL':
                    edit_list.append(value)
            # Check if the number is possible
            if len(edit_list) < number:
                number = len(edit_list)
            return edit_list, number

    def write_file(self, name_file, list_args):
        """
        Writes the arguments (list_args) to a file.
        :param name_file:   File name (and path) after which the arguments are written
        :param list_args:   list of arguments with the following arguments:
                            mkdir_path:     Path folder to be created
                            arg_mutect2:    Argument to run Mutect2
                            file_output:    Path and first part of output files

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


def make_objects(df_selection):
    """
    Creates the dictionary with the objects per participant/sample
    :param df_selection:    the custom data frame from which to create objects and a dictionary
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


def arguments_to_file(dict_samples, head_path, chrom, number_of_tumors=None, number_of_hc=None, type_sample='both',
                      type_aln='bowtie', type_aln2='bowtie2'):
    """
    Ensures that all arguments are created and written.
    :param dict_samples:     a dictionary containing the sample_num as key and as value Sample objects
    :param head_path:       Main path to where the file will be saved
    :param number_of_tumors: Number of hc you want to combine while running Mutect2
    :param number_of_hc:    Number of hc you want to combine while running Mutect2
    :param type_sample:     What type of tumor you want to have ("both", "tFL" or "FL")
    :param type_aln:
    :param type_aln2:
    :return:
    """
    # Loop over keys from dict_sample
    compare_hc_tum = open(f'{head_path}{chrom}_compare_hc_tumor_{type_aln}_{type_sample}.txt', "w")
    manual_comparison = open(f'{head_path}{chrom}_manual_comparison_{type_aln}_{type_sample}.txt', "w")
    mutect2_comparison = open(f'{head_path}{chrom}_mutect2_comparison_{type_aln}_{type_sample}.txt', "w")
    for key in dict_samples:
        dict_samples[key].get_arguments(head_path, chrom, compare_hc_tum, manual_comparison, mutect2_comparison,
                                        number_of_tumors, number_of_hc,
                                        type_sample, type_aln,
                                        type_aln2)
    compare_hc_tum.close()
    manual_comparison.close()
    mutect2_comparison.close()


def main():
    """

    :return:
    """
    # the path to the file
    path_file = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/EGAD00001000292_metadata/delimited_maps/Sample_File.map"
    # Main path to where the file will be saved
    head_path = sys.argv[1]
    # Number of tumors you want to combine while running Mutect2.
    number_of_tumors = int(sys.argv[2])
    # Number of hc you want to combine while running Mutect2
    number_of_hc = int(sys.argv[3])
    # What type of tumor you want to have ("both", "tFL" or "FL")
    type_sample = sys.argv[4]
    # Which method was used
    # Name of the folder
    type_aln = sys.argv[5]  # bowtie, bwa_aln, bwa_mem
    # Piece with which the file name begins
    type_aln2 = sys.argv[6]  # bowtie2, aln, mem
    chrom = sys.argv[7]
    # Filter file
    df_selection = filter_file(path_file)
    # Make objects for all samples/participants
    dict_samples = make_objects(df_selection)

    arguments_to_file(dict_samples, head_path, chrom, number_of_tumors, number_of_hc, type_sample, type_aln, type_aln2)


if __name__ == "__main__":
    main()
