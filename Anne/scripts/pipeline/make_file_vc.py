#!/usr/bin/env python3
import pandas as pd
from itertools import combinations
import sys

from Sample import Sample
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config




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
    compare_hc_tum = open(f'{head_path}{chrom}_compare_hc_tumor_{type_aln}_{type_sample}.txt', "w")
    manual_comparison = open(f'{head_path}{chrom}_manual_comparison_{type_aln}_{type_sample}.txt', "w")
    mutect2_comparison = open(f'{head_path}{chrom}_mutect2_comparison_{type_aln}_{type_sample}.txt', "w")

    # Loop over keys from dict_sample
    for key in dict_samples:
        # dict_samples[key] > Sample()
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
    config = get_config()
    # the path to the file
    path_file = config['sample_file']
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
