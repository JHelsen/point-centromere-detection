# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 11:08:59 2025

@author: Kausthubh R
"""

import os
from tqdm import tqdm
from Bio import SeqIO
from copy import deepcopy

# =============================================================================
# Please make sure to run this script in the directory created in the final
# step of 4_reverse_blast_analysis_for_species_tree.py script i.e. the 
# directory whose name ends with "_seqs_for_tree_building"
# =============================================================================

base_list = os.listdir()
base_path = os.getcwd()


#%% getting alignments for tree building and concatenation

alignments = {}
for i in tqdm(base_list, position=0, leave=True):
    os.chdir(i)
    try:
        alignments[i] = [y for y in SeqIO.parse([x for x in os.listdir() if x.endswith("trim.phy")][0], "fasta")]
    except IndexError or "NewickError":
        print(f"\n{i} - no alignment\n")
    os.chdir(base_path)
del i

#%% renaming alignment sequences
# =============================================================================
# This cell renames each sequence in the alignment. Each sequence name is 
# replaced with the corresponding species name
# =============================================================================
alignments_renamed = {}
for i in tqdm(alignments, position=0, leave=True):
    temp = deepcopy(alignments[i])
    for j in temp:
        name = j.id[:j.id.rfind(".")]
        name = name[name.find("_")+1:]
        j.id = name
        j.description = ""
    alignments_renamed[i] = temp
del i, temp, j, name
        

#%% writing the alignments out
path2alignments = input("Please enter the path to the empty directory where you want to write the modified alignments out:\n")
os.chdir(path2alignments)

for i in tqdm(alignments_renamed, leave=True, position=0):
    with open(f"{i}_renamed_align_trim.phy", 'w') as f:
        [SeqIO.write(x, f, 'fasta') for x in alignments_renamed[i]]
        f.close()
del i, f
