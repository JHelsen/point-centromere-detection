# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 13:21:32 2024

@author: Kausthubh R
"""
import pandas as pd
from Bio import SeqIO
import copy
import os
from tqdm import tqdm

#obj1 is the name of the tblout file
def read_phmmer_tblout(obj1):
    temp = open(obj1, 'r')
    temp_con = temp.read().split("\n")[0:-10]
    temp1 = []
    for i in range(len(temp_con[1:])):
        temp2 = copy.deepcopy(temp_con[i+1])
        temp2 = temp2.split(" ")
        temp2 = list(filter(None, temp2))
        temp1.append(temp2)
    del i, temp2
    for i in temp1:
        try:
            for j in range(0,18):
                try:
                    i[j] = float(i[j])
                except ValueError:
                    continue
        except IndexError:
            continue
    del i, j, temp1[1], temp1[0][0]
    #%
    temp1[0][3] = temp1[0][3] + " " + temp1[0][4]
    temp1[0][0] = temp1[0][0] + " " + temp1[0][1]
    temp1[0][20] = temp1[0][20] + " " + temp1[0][21] + " " + temp1[0][22]
    temp1[0][6] = "full sequence " + temp1[0][6]
    temp1[0][7] = "full sequence " + temp1[0][7]
    temp1[0][8] = "full sequence " + temp1[0][8]
    temp1[0][9] = "best 1 domain " + temp1[0][9]
    temp1[0][10] = "best 1 domain " + temp1[0][10]
    temp1[0][11] = "best 1 domain " + temp1[0][11]
    temp1[0][12] = "domain number estimation " + temp1[0][12]
    temp1[0][13] = "domain number estimation " + temp1[0][13]
    temp1[0][14] = "domain number estimation " + temp1[0][14]
    temp1[0][15] = "domain number estimation " + temp1[0][15]
    temp1[0][16] = "domain number estimation " + temp1[0][16]
    temp1[0][17] = "domain number estimation " + temp1[0][17]
    temp1[0][18] = "domain number estimation " + temp1[0][18]
    temp1[0][19] = "domain number estimation " + temp1[0][19]
    
    del  temp1[0][22], temp1[0][21], temp1[0][4], temp1[0][1], temp1[int(len(temp1)-1)]
    
    for i in temp1[1:]:
        for j in range(19, len(i)):
            i[18] = i[18] + " " + i[j]
            i[j] = []
        
    del i, j
    for i in range(len(temp1)):
        temp1[i] = [ele for ele in temp1[i] if ele != []]
    del i
    
    data = pd.DataFrame(temp1[1:len(temp1)], columns=temp1[0])
    return data

cerevisiae_proteome = {x.id : x for x in SeqIO.parse(r"Z:\Kausthubh\1_postdocs\Jana\Centromere_project_phylogeny\cerevisiae_proteome\Sc_Proteome.fasta", 'fasta')}

tblout_list = os.listdir()
#%%
# =============================================================================
# Importing best (1 domain) scoring hit from each phmmer search
# =============================================================================
print("Importing best (1 domain) scoring hit from each phmmer search\n")
selected_seqs = {}
repeats = {}
for i in tqdm(tblout_list):
    temp = read_phmmer_tblout(i)
    selected_id = temp[temp["best 1 domain score"] == max(list(temp["best 1 domain score"]))]["target name"]
    if len(selected_id) > 1:
        input("wait")    
    selected_seqs[i] = cerevisiae_proteome[selected_id.iloc[0]]
del temp, selected_id, i

# =============================================================================
# Identifying those phmmer searches in which the best (1 domain) scoring hit
# is identified in some other phmmer search
# =============================================================================

print("Identifying duplicate hits\n")
for i in tqdm(selected_seqs):
    id_check = [selected_seqs[x].id for x in selected_seqs if x != i]
    if selected_seqs[i].id in id_check:
        repeats[i] = selected_seqs[i].id
del i, id_check

repeats = dict(sorted(repeats.items(), key=lambda item: item[1]))

# =============================================================================
# Clustering the identified hits
# =============================================================================
repeat_clusters = {}
completed = []
for i in repeats.values():
    if i in completed:
        continue
    repeat_clusters[i] = [x for x in repeats if repeats[x] == i]
    completed.append(i)
del i

# =============================================================================
# Selecting the first element in each cluster
# =============================================================================
repeats_selects = [x[0] for x in repeat_clusters.values()]
repeats_rejects = [x for x in repeats if x not in repeats_selects]

# =============================================================================
# Writing a log file documenting which markers had identical best (1 domain)
# hits, and which markers were not included for further analysis
# =============================================================================
print("Writing a log file documenting duplicates\n")
with open(os.path.join(os.path.dirname(os.getcwd()), "duplicate_markers.txt"), 'w') as f:
    [f.write(f"{x}\t{repeats[x]}\n") for x in repeats]
    f.write("\n")
    f.write("The following markers were not used since they returned a best (1 domain) scoring hit identical to one of the previously selected markers\n")
    [f.write(f"{x}\n") for x in repeats_rejects]
    # f.write(f"{datetime.now()} - Current,Peak\n")
    f.close()
    
# =============================================================================
# Writing the hits out into separate fasta files in a separate folder
# =============================================================================
print("Writing selected sequences into a separate folder\n")
os.chdir(os.path.dirname(os.getcwd()))
os.mkdir(f"marker_seqs_{len([x for x in selected_seqs if x not in repeats_rejects])}")
os.chdir(f"marker_seqs_{len([x for x in selected_seqs if x not in repeats_rejects])}")

for i in tqdm(selected_seqs):
    if i in repeats_rejects:
        continue
    # input()
    SeqIO.write(selected_seqs[i], f"{selected_seqs[i].id}.fasta", 'fasta')