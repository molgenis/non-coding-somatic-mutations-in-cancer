"""
Checken of er normal tissue in zit
"""

import gzip

specimen = gzip.open("D:/Hanze_Groningen/STAGE/NEW PLAN/all_specimen.tsv.gz", 'rt')
header = specimen.readline().strip().split("\t")
print(header)
samples = {}

for line in specimen:
    line = line.strip()
    elems = line.split("\t")
    sample_id = elems[0]
    annotation = elems[6]
    cancer = 1
    if "Normal" in annotation:
        cancer = 0
    samples[sample_id] = cancer

specimen.close()

project_file = gzip.open("D:/Hanze_Groningen/STAGE/NEW PLAN/ALL-US.tsv.gz", 'rt')
header = project_file.readline().strip().split("\t")

ctr = 0
sampleIdsSeen = set()

for line in project_file:
    line = line.strip()
    elems = line.split("\t")
    sample_id = elems[3]
    is_cancer = samples.get(sample_id)

    if elems[33] == 'WGS':
        if is_cancer is not None and sample_id not in sampleIdsSeen:
            if is_cancer == 1:
                # this is a cancer sample.
                ctr += 1
            else:
                print(sample_id)
                # this is not a cancer sample
        else:
            x = 'x'
            # if elems[33] != 'WGS':
            #     print(elems[33])
            #     print(f'---- {sampleid}')
            # specimen ID not found in annotation file
        sampleIdsSeen.add(sample_id)

project_file.close()

print("{} ids were cancer id out {}".format(ctr, len(sampleIdsSeen)))
