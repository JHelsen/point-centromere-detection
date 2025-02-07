# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 19:33:42 2024

@author: Kausthubh R
"""

import random
from Bio import SeqIO
import os

# =============================================================================
# Loading alignment as a fasta file
# =============================================================================
path2supermatrixalignment = input("Please provide the path to the supermatrix alignment. Please include the alignment filename and extension:\n")
alignment = [x for x in SeqIO.parse(path2supermatrixalignment, 'fasta')]

# =============================================================================
# Deciding subsample size
# =============================================================================
subsample = 10000
range_1 = range(int(sum([len(x) for x in alignment])/len(alignment)))

# =============================================================================
# Subsampling using the random module
# =============================================================================
subsample_1 = random.sample(range_1, subsample)
subsample_2 = random.sample(range_1, subsample)
subsample_3 = random.sample(range_1, subsample)

#%%
# =============================================================================
# Checking intersections between the different subsamples
# =============================================================================
inter_12 = set(subsample_1) & set(subsample_2)
inter_23 = set(subsample_3) & set(subsample_2)
inter_31 = set(subsample_1) & set(subsample_3)

#%%

# =============================================================================
# Extracting subsample of the sequences based on the subampled range
# =============================================================================
subsample_1_sequence = {}
subsample_2_sequence = {}
subsample_3_sequence = {}
for i in alignment:
    subsample_1_sequence[i.id] = "".join([i.seq[x] for x in subsample_1])
    subsample_2_sequence[i.id] = "".join([i.seq[x] for x in subsample_2])
    subsample_3_sequence[i.id] = "".join([i.seq[x] for x in subsample_3])
    
#%%
# =============================================================================
# Writing the sequences out into fasta files
# =============================================================================

path2subsample_directory = input("Please provide the path to the directory where the subsample files need to be written:\n")
os.chdir(path2subsample_directory)

with open("subsample_1_sequence.fasta", 'w') as f:
    [f.write(f">{x}\n{subsample_1_sequence[x]}\n") for x in subsample_1_sequence]
    f.close()
    
with open("sites_subsampled_1.txt", 'w') as f:
    f.write("Subsampled based on Python indexing\n")
    [f.write(f"{x},") for x in subsample_1]
    f.close()

with open("subsample_2_sequence.fasta", 'w') as f:
    [f.write(f">{x}\n{subsample_2_sequence[x]}\n") for x in subsample_2_sequence]
    f.close()
    
with open("sites_subsampled_2.txt", 'w') as f:
    f.write("Subsampled based on Python indexing\n")
    [f.write(f"{x},") for x in subsample_2]
    f.close()
    
with open("subsample_3_sequence.fasta", 'w') as f:
    [f.write(f">{x}\n{subsample_3_sequence[x]}\n") for x in subsample_3_sequence]
    f.close()
    
with open("sites_subsampled_3.txt", 'w') as f:
    f.write("Subsampled based on Python indexing\n")
    [f.write(f"{x},") for x in subsample_3]
    f.close()