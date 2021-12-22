import pandas as pd
import glob

path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/'
path_files = f'{path}down*.tsv'
# Loop over all files
set_donor = set()
list_donor = []
for fname in glob.glob(path_files):
    print(fname)
    # Read file
    df = pd.read_csv(fname, sep='\t')
    set_donor.update(set(df['icgc_donor_id']))
    list_donor.extend(list(set(df['icgc_donor_id'])))

print(len(set_donor))
print(len(list_donor))
