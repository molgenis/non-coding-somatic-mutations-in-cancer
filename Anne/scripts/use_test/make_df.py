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
# Select the columns from the column named 'FORMAT'
df_only_samples = df.iloc[:, df.columns.get_loc("FORMAT"):]
print(df_only_samples.columns)

# Create a set of all abbreviations in the FORMAT column that are separated by :
set_format = set()
# Loop over each line of the vcf file
for index, row in df_only_samples.iterrows():
    list_format = row['FORMAT'].split(':')
    set_format.update(list_format)
print(set_format)

# Make a dataframe with column names equal to the values from the set_format
df_format_value = pd.DataFrame(columns=list(set_format))
# Loop over each line of the vcf file
for index, row in df_only_samples.iterrows():
    # Create list that will contain values for each row
    list_values = list()
    # Split column with name FORMAT on :
    names_format = row['FORMAT'].split(':')
    # Split the last (also the second column) on :
    values_format = row[1].split(':')
    # Loop over the column names of the new dataframe (set_format)
    for col_name in list(set_format):
        # Check if the col_name occurs in the abbreviations of the FORMAT column (in this specific row)
        # It is possible that one row contains some values and the other row does not
        if col_name in names_format:
            # Find the index the specific abbreviation (col_name) in names_format
            index_value = names_format.index(col_name)
            # Add the specific value corresponding to the index (index_value) in the values_format list
            list_values.append(values_format[index_value])
        else:
            # If the col_names does not exist in the names_format, it means that this row
            # does not have that particular value. NaN is then added to list_values
            list_values.append(float("NaN"))
    # Always add the list_values as the last row in the new dataframe
    df_format_value.loc[len(df_format_value)] = list_values
# Save the dataframe as a csv file
df_format_value.to_csv(f'{sys.argv[2]}{df_only_samples.columns[-1].split(".")[0]}.csv', sep='\t', encoding='utf-8')
