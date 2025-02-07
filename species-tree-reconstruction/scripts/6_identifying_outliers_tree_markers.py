# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 11:39:36 2025

@author: Kausthubh R
"""


import ete3
import os
from tqdm import tqdm
from statistics import median
import shutil
from Bio import SeqIO

base_list = os.listdir()
base_path = os.getcwd()

# importing trees
trees = {}
for i in tqdm(base_list, position=0, leave=True):
    os.chdir(i)
    try:
        trees[i] = ete3.Tree([x for x in os.listdir() if x.endswith(".tree")][0], format=1)
    except IndexError or "NewickError":
        print(f"\n{i} - no tree\n")
    os.chdir(base_path)
del i

#%% checking outliers
# =============================================================================
# Outliers are detected if their branch length exceeds 20 times the median
# branch length. If the median branch length is less than 1e-8, no action is 
# taken and the user is asked to check it manually.
# 
# If more than 33# of the leaves are classified as potential outliers, no
# action is taken and the user is asked to check it manually.
# =============================================================================
outgroups = ["Wickerhamomyces_anomalus","Yarrowia_lipolytica","Pichia_kudriavzevii","Candida_albicans"]
unusual = []
outliers = {}
for i in tqdm(trees, position=0, leave=True):
    
    nodes = {x.name : [x.dist, [y for y in x.get_children() if y.is_leaf()]] for x in trees[i].traverse("preorder") if not [y for y in outgroups if y in x.name]}
    median_point = median([x[0] for x in nodes.values()])
    if median_point < 1e-8:
        print(f"\nThe median is too small (<1e-8) in {i}. Please check manually\n")
        continue
   
    potential_outliers = [x for x in nodes if nodes[x][0] > 20*median_point]
    if len(potential_outliers) > len(nodes)/3:
        print(f"\nUnusually high number of outliers (more than 33% of existing leaves) in {i}. Please check manually\n")
        unusual.append(i)
        continue
    if potential_outliers:
        outliers[i] = potential_outliers
    
del i, median_point, potential_outliers, nodes


for i in outliers:
    if i in unusual:
        continue
    for j in outliers[i]:
        try:
            float(j)
            unusual.append(i)
            
            break
        except ValueError:
            unusual.append(i)
            continue
del i, j

outliers = {x : outliers[x] for x in outliers if x not in unusual}
#%% removing outliers and preparing for fresh tree building
outliers_final = []
if outliers:
    for i in tqdm(outliers, position=0, leave=True):
        os.chdir(i)
        count = len([x for x in os.listdir() if x.startswith("old")])
        
        seqs = [x for x in SeqIO.parse([y for y in os.listdir() if y.endswith(".fasta") if "small" not in y and "updated" not in y][0], 'fasta')]
        
        outlier_seqs = [x for x in seqs if x.id in outliers[i]]
        seqs = [x for x in seqs if x.id not in outliers[i]]
        
        with open(f"{i}_outliers_{count}.fasta", 'w') as f:
            [SeqIO.write(x, f, 'fasta') for x in outlier_seqs]
            f.close()
        files = [x for x in os.listdir() if not x.startswith("old")]
        
        os.mkdir(f"old_{count}")
        dest_path = os.path.join(os.getcwd(), f"old_{count}")
        [shutil.move(x, os.path.join(dest_path, x)) for x in files]
        
        with open(f"{i}_homologs.fasta", 'w') as g:
            [SeqIO.write(x, g, 'fasta') for x in seqs]
            g.close()
        os.chdir(base_path)
        outliers_final.append(i)
        # input("hi")

if unusual:
    for i in tqdm(unusual, position=0, leave=True):
        os.chdir(i)
        count = len([x for x in os.listdir() if x.startswith("old")])
        
        try:
            outlier_seqs = [x.name for x in ete3.Tree([x for x in os.listdir() if "removed" in x][0], format=1) if x.is_leaf()]
        except IndexError:
            print(f"\nNo outliers found in {i}\n")
            os.chdir(base_path)
            continue
        outliers_final.append(i)
        seqs = [x for x in SeqIO.parse([y for y in os.listdir() if y.endswith(".fasta") if "small" not in y and "updated" not in y][0], 'fasta')]
        outlier_seqs = [x for x in seqs if x.id in outlier_seqs]        
        seqs = [x for x in seqs if x.id not in [y.id for y in outlier_seqs]]
        
        with open(f"{i}_outliers_{count}.fasta", 'w') as f:
            [SeqIO.write(x, f, 'fasta') for x in outlier_seqs]
            f.close()
        files = [x for x in os.listdir() if not x.startswith("old")]
        
        os.mkdir(f"old_{count}")
        dest_path = os.path.join(os.getcwd(), f"old_{count}")
        [shutil.move(x, os.path.join(dest_path, x)) for x in files]
        
        with open(f"{i}_homologs.fasta", 'w') as g:
            [SeqIO.write(x, g, 'fasta') for x in seqs]
            g.close()
        os.chdir(base_path)

else:
    print("\nNo outliers found\n")

for i in outliers_final:
    print(i)
    