# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:54:31 2024

@author: Kausthubh R
"""

import os
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from joblib import Parallel, delayed

# =============================================================================
# It is recommended to commence your run in the same directory which contains 
# the results of the 1_parallelized_tblastn_searches.py script, the 
# 2_tblastn_results_analysis_orf_extraction.py script, and the 
# 3_reverse_blastp_against_cerevisiae.py script
# =============================================================================

path2genomes = input("Please enter the path to the directory containing your genomes for the analysis:\n")

listofgenomes = os.listdir(path2genomes)

run_title = os.path.basename(os.getcwd())

os.chdir(f"{run_title}_selected_orfs_reciprocal_blast")

base_list = [x for x in os.listdir() if not x.startswith(".")]
base_path = os.getcwd() 

#%% importing reciprocal blastp results
print("Importing reciprocal BLASTp results\n")
base_list_import = [base_list[i:i + int(os.cpu_count())] for i in range(0, len(base_list), int(os.cpu_count()))]
def import_blastp(i:str, listofgenomes:list, base_path:str):
    os.chdir(i)
    orfs = {}
    # print(f"Importing reciprocal blastp results for {i}\n")
    temp2 = pd.DataFrame(columns=["query id","subject id", "evalue"], index = listofgenomes)
    
    for j in os.listdir():
        name = j[j.find("_")+1:j.rfind("reciprocal")-1] + ".fna" #only for marker seqs for species tree
        
        file = open(j, 'r').read().split("\n")
        # number_of_processed_queries = int(file[-2][file[-2].rfind("queries")-2])
        inds = [x+1 for x in range(len(file)) if "hits found" in file[x]]
        query_inds = [x for x in range(len(file)) if "Query" in file[x]]
        temp1 = []
        temp3 = []
        temp4 = []
        for k in range(len(inds)):
            if file[inds[k]].split("\t")[0] == "# BLASTP 2.5.0+" or len(file[inds[k]].split("\t")) == 1:
                temp1.append("No reciprocal hit found")
                temp3.append("No reciprocal hit found")
                temp4.append(file[query_inds[k]][file[query_inds[k]].find(":")+2:])
                continue
            temp1.append(float(file[inds[k]].split("\t")[0]))
            temp3.append(file[inds[k]].split("\t")[-2])
            temp4.append(file[query_inds[k]][file[query_inds[k]].find(":")+2:])
        
        # =============================================================================
        # Only imports those hits which have mapped back to the corresponding starting 
        # ID in S. cerevisiae. This ONLY works for marker sequences used for species
        # tree reconstruction i.e. when the folder names match the SGD protein ID        
        # =============================================================================
        temp_zipped = list(zip(temp4, temp3, temp1))
        temp_zipped = [z for z in temp_zipped if z[1] == i]
        
        if len(temp_zipped) > 1:
            temp2["evalue"][name] = [x[2] for x in temp_zipped]
            temp2["query id"][name] = [x[0] for x in temp_zipped]
            temp2["subject id"][name] = [x[1] for x in temp_zipped]
        if len(temp_zipped) == 1:
            temp2["evalue"][name] = temp_zipped[0][2]  
            temp2["query id"][name] = temp_zipped[0][0]
            temp2["subject id"][name] = temp_zipped[0][1]
            
        
    orfs[i] = temp2.fillna("No hit found")
    os.chdir(base_path)
    return orfs

orfs_by_folder = []
for i in tqdm(base_list_import, position=0, leave=True):
    temp = Parallel(n_jobs=int(os.cpu_count()))(delayed(import_blastp)(j, listofgenomes, base_path) for j in i)
    [orfs_by_folder.append(x) for x in temp]
del i, temp, base_list_import

orfs = {k: v for d in orfs_by_folder for k, v in d.items()}
print("Complete\n")

#%% organizing results as a dataframe and checking if there are any markers absent in more than 50% of species
print("Organizing imported results into a dataframe\n")
verified_profile = pd.DataFrame("No hit found",index=listofgenomes, columns=orfs.keys())
verified_profile_numbers = pd.DataFrame(0,index=listofgenomes, columns=orfs.keys())
for i in orfs.keys():
    for j in listofgenomes:
        if type(orfs[i]["query id"][j]) is list:
            verified_profile[i][j] = orfs[i]["query id"][j]
            verified_profile_numbers[i][j] = len(orfs[i]["query id"][j])
            continue
        if orfs[i]["query id"][j] == "No hit found":
            continue
        else:
            verified_profile[i][j] = orfs[i]["query id"][j]
            verified_profile_numbers[i][j] = 1
del i, j
#% checking mostly missing profiles

print("Dropping those sequences which do not have homologs in at least 50% of the species\n")
check = [x for x in verified_profile_numbers if list(verified_profile_numbers[x]).count(0) >= len(listofgenomes)/2]

[verified_profile.drop(x, axis=1, inplace=True) for x in check]
[verified_profile_numbers.drop(x, axis=1, inplace=True) for x in check]
print("Complete\n")
#%% importing the protein orf predicted sequences
print("Importing the predicted ORF sequences\n")
#windows
os.chdir(os.path.dirname(os.getcwd()))
os.chdir(f"{run_title}_selected_orfs")

orf_seqs = []

def import_orf_seqs(i:str, verified_profile:pd.DataFrame()):
    os.chdir(i)
    
    compare_list = list(verified_profile[i])
    compare_list = [item for sublist in compare_list for item in (sublist if isinstance(sublist, list) else [sublist])]
    compare_list = [x for x in compare_list if x != "No hit found"]
    
    temp = [[x for x in SeqIO.parse(y, 'fasta') if [z for z in compare_list if x.id in z]] for y in os.listdir()]
    temp = [x for y in temp for x in y]
    
    x = {}
    
    x[i] = temp
    os.chdir(os.path.dirname(os.getcwd()))

    return x

temp = [x for x in base_list if x not in check]    
base_list_import = [temp[i:i + int(os.cpu_count())] for i in range(0, len(temp), int(os.cpu_count()))]
for i in tqdm(base_list_import, position=0, leave=True):
    temp = Parallel(n_jobs=int(os.cpu_count()))(delayed(import_orf_seqs)(j, verified_profile) for j in i)
    [orf_seqs.append(x) for x in temp]
del i, temp, base_list_import

orf_seqs = {k: v for d in orf_seqs for k, v in d.items()}
print("Complete\n")
#%% selecting at most one sequence per species regardless of copy number
orf_seqs_to_write = {}
for i in tqdm(orfs, position=0, leave=True):
    if i in check:
        continue
    to_write = []
    for j in listofgenomes:
        temp = [x for x in orf_seqs[i] if j[:j.rfind(".")] in x.id]
        if not temp:
            continue
        if len(temp) == 1:
            to_write.append(temp[0])
        if len(temp) > 1:
            temp = [x for x in temp if len(x) == max([len(y) for y in temp])]
            to_write.append(temp[0])
    orf_seqs_to_write[i] = to_write
del temp, i, to_write, j

#%% writing selected orfs out
print("Writing selected, reciprocally verified ORFs out\n")
os.chdir(os.path.dirname(os.getcwd()))
verified_profile.to_csv(f"{run_title}_reciprocally_verified_orfs.csv")
verified_profile_numbers.to_csv(f"{run_title}_reciprocally_verified_orfs_counts.csv")

os.mkdir(f"{run_title}_seqs_for_tree_building")
os.chdir(f"{run_title}_seqs_for_tree_building")

for i in tqdm(orf_seqs_to_write, position=0, leave=True):
    os.mkdir(i)
    os.chdir(i)
    with open(f"{i}_homologs.fasta", 'w') as f:
        [SeqIO.write(x, f, 'fasta') for x in orf_seqs_to_write[i]]
        f.close()
    os.chdir(os.path.dirname(os.getcwd()))
del i
print("Complete/n")