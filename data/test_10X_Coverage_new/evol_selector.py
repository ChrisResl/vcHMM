import os
import pysam


location = "/home/alexander/Documents/uni/sequence/vi_HMM/vcHMM/data/test_10X_Coverage_new"
goodsam = []

for file in os.listdir(location):
    if file.endswith('.sam'):
        samfile = pysam.AlignmentFile(file, 'rb')
        sam = [[read.reference_start, read.query_name, read.query_sequence, read.cigartuples,
                read.mapping_quality, read.query_qualities.tolist()] for read in samfile]
        if any(read[3] is None for read in sam):
            continue
        else:
            goodsam.append(os.path.join(file))

goodsam.sort()

for file in goodsam:
    print(file)
