#!/usr/bin/env python3
import pandas as pd
import io
import sys


# https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
def read_vcf(path):
    """

    :param path:
    :return:
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


# Read vcf file
df = read_vcf(sys.argv[1].strip())
# Make list of the last two columns. The first index is always hc and the second index is always tumor sample
hc_tum = df.columns[-2:]
# Select the columns from the column named 'FORMAT'
df_only_samples = df.iloc[:, df.columns.get_loc("FORMAT"):]

# Create a set of all abbreviations in the FORMAT column that are separated by :
set_format = set()
# Loop over each line of the vcf file
for index, row in df_only_samples.iterrows():
    list_format = row['FORMAT'].split(':')
    set_format.update(list_format)
print(set_format)

# Make a dataframe with column names equal to the values from the set_format
df_format_value_hc = pd.DataFrame(columns=list(set_format))
df_format_value_tum = pd.DataFrame(columns=list(set_format))
# Loop over each line of the vcf file
for index, row in df_only_samples.iterrows():
    # Create list that will contain values for each row
    list_values_hc = list()
    list_values_tum = list()
    # Split column with name FORMAT on :
    names_format = row['FORMAT'].split(':')
    # Split the hc and tum column on :
    values_hc_format = row[hc_tum[0]].split(':')
    values_tum_format = row[hc_tum[1]].split(':')
    # Loop over the column names of the new dataframe (set_format)
    for i in list(set_format):
        # Check if the col_name occurs in the abbreviations of the FORMAT column (in this specific row)
        # It is possible that one row contains some values and the other row does not
        if i in names_format:
            # Find the index the specific abbreviation (col_name) in names_format
            index_value = names_format.index(i)
            # Add the specific value corresponding to the index (index_value) in the values_format list
            list_values_hc.append(values_hc_format[index_value])
            list_values_tum.append(values_tum_format[index_value])
        else:
            # If the col_names does not exist in the names_format, it means that this row
            # does not have that particular value. NaN is then added to list_values
            list_values_hc.append(float("NaN"))
            list_values_tum.append(float("NaN"))
    # Always add the list_values as the last row in the new dataframe
    df_format_value_hc.loc[len(df_format_value_hc)] = list_values_hc
    df_format_value_tum.loc[len(df_format_value_tum)] = list_values_tum
# Save the dataframe as a csv file
df_format_value_hc.to_csv(f'{sys.argv[2]}{hc_tum[0].split(".")[0]}_{hc_tum[1].split(".")[0].split("_")[1]}_hc.csv',
                          sep='\t', encoding='utf-8')
df_format_value_tum.to_csv(f'{sys.argv[2]}{hc_tum[1].split(".")[0]}.csv', sep='\t', encoding='utf-8')
