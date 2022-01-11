import pandas as pd
import glob

path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/'
path_files = f'{path}down*.tsv'
# Loop over all files
set_donor = set()
list_donor = []
set_samples = set()
set_matched_sample = set()
list_samples = []
list_matched_sample = []

for fname in glob.glob(path_files):
    print(fname)
    # Read file
    df = pd.read_csv(fname, sep='\t')
    print(df.columns)
    set_donor.update(set(df['icgc_donor_id']))
    list_donor.extend(list(set(df['icgc_donor_id'])))
    print(len(set_donor))
    print(len(list_donor))
    print('-----')
    set_samples.update(set(df['icgc_sample_id']))
    list_samples.extend(list(set(df['icgc_sample_id'])))
    print(len(set_donor))
    print(len(list_donor))
    print('-----')
    set_matched_sample.update(set(df['matched_icgc_sample_id']))
    list_matched_sample.extend(list(set(df['matched_icgc_sample_id'])))
    print(len(set_donor))
    print(len(list_donor))

print('END')
print(len(set_donor))
print(len(list_donor))
print('-----')
print(len(set_samples))
print(len(list_samples))
print('-----')
print(len(set_matched_sample))
print(len(list_matched_sample))
print('ENDENDENDENDEND')
