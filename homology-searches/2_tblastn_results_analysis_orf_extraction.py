# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 13:31:13 2023

@author: Kausthubh R
"""

import os
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import BiopythonWarning
import pandas as pd
# from statistics import mean
from orffinder import orffinder
# from math import isnan
from tqdm import tqdm
from joblib import Parallel, delayed
from itertools import islice
import tracemalloc
import warnings
warnings.filterwarnings("ignore")



if os.name == 'nt':
    def winpath2linuxpath(obj1):
        temp = obj1
        temp = temp.replace("\\", "/")
        temp = "/mnt/" + temp[0].lower() + "/" + temp[3:]
        return temp
    
    def serverpath2winpath(obj1):
        temp = obj1
        temp = r"Z:" + temp[6:]
        temp = temp.replace("/", r"\\")
        return temp
    
    def winpath2serverpath(obj1):
        temp = obj1
        temp = temp.replace("Z:", "/g/dey")
        temp = temp.replace("\\", "/")
        temp = temp.replace("//", "/")


# =============================================================================
# This block of code reads tblastn search results [command used - tblastn -db 
# {list_of_blast_dbs[j]} -query {i} -outfmt '7 evalue bitscore score length 
# pident nident qstart qend sstart send sseqid sseq' -out {name}_{j}.txt] and 
# makes a dictionary of type Protein:Dataframe
# =============================================================================
path2genomes = input("Please enter the path to the directory containing your genomes for the analysis:\n")

folder_name = input("Please enter the path to the tblastn search results:\n")

list_of_genomes = dict([(x, os.path.join(path2genomes, x)) for x in os.listdir(path2genomes)])
print("\nLoading genomes\n")
genomes = {}
for i in tqdm(list_of_genomes, unit=" genomes", position=0, leave=True):
    genomes[i] = [x for x in SeqIO.parse(list_of_genomes[i], 'fasta')]
    # print(i)
del i

# genomes = list(genomes.keys())
print("\nComplete\n")
#%% importing tblastn results
# =============================================================================
# This block of code only travels to each directory, opens the tblastn results
# in a dataframe, and only retains those hits whose evalue is lower than 1e-10.
# In case there are no hits with an evalue < frst_evalue_thresold, only the hit
# with the lowest evalue is selected
# =============================================================================

print("\nImporting tblastn results\n")
os.chdir(folder_name)
frst_evalue_thresold = 1e-10
selected_scaffolds = {}
base_list = [x for x in os.listdir() if not x.startswith(".")]
base_list_import = [base_list[i:i + int(os.cpu_count())] for i in range(0, len(base_list), int(os.cpu_count()))]
base_path = os.getcwd() #the starting folder should be the folder containing the tblastn search results
def get_best_tblastn_result(i, base_path, frst_evalue_thresold):
    # i = base_list[i]
    os.chdir(i)
    temp_dict= {}
    # print(f"Importing scaffolds for {i}\n")
    temp2 = pd.DataFrame(columns=["subject id", "s. start", "s. end", "evalue"], index = os.listdir())
    # temp2 = []
    for j in os.listdir():
       
        file = open(j, 'r').read().split("\n")
        if len(file) <= 6:
            temp2["subject id"][j] = "No hit found"            
            continue
        
        # seq_name = j[:j.find("YEAST_")] + "YEAST"
                
        headers = file[3]
        headers = headers[headers.find(":")+1:].split(", ")
        headers[0] = "evalue"
        
        values = file[5:len(file)-2]
        values = [x.split("\t") for x in values]
        values = [list(map(float, [x[y] for y in range(len(x)) if not y > 9])) + [x[y] for y in range(len(x)) if y > 9]  for x in values]
        results = pd.DataFrame(values, columns = headers)
        # input()
        temp = results[results["evalue"] < frst_evalue_thresold]

        if len(temp) == 0:
            temp_2 = results[results["evalue"] == min(list(results["evalue"]))]
            temp2["subject id"][j] = temp_2["subject id"][temp_2.index[0]]
            temp2["s. start"][j] = int(temp_2["s. start"][temp_2.index[0]])
            temp2["s. end"][j] = int(temp_2["s. end"][temp_2.index[0]])
            temp2["evalue"][j] = float(temp_2["evalue"][temp_2.index[0]])
            continue
        
        if len(temp) == 1:
            temp_2 = temp[temp["evalue"] == min(list(temp["evalue"]))]
            temp2["subject id"][j] = temp_2["subject id"][temp_2.index[0]]
            temp2["s. start"][j] = int(temp_2["s. start"][temp_2.index[0]])
            temp2["s. end"][j] = int(temp_2["s. end"][temp_2.index[0]])
            temp2["evalue"][j] = float(temp_2["evalue"][temp_2.index[0]])
            continue
        
        temp2["subject id"][j] = [temp["subject id"][a] for a in temp.index]
        temp2["s. start"][j] = [int(temp["s. start"][a]) for a in temp.index]
        temp2["s. end"][j] = [int(temp["s. end"][a]) for a in temp.index]
        temp2["evalue"][j] = [float(temp["evalue"][a]) for a in temp.index]
            
    temp_dict[i] = temp2
    os.chdir(base_path)
    return temp_dict

selected_scaffolds = []
for i in tqdm(base_list_import, unit=" sequences", position=0, leave=True):
    temp = Parallel(n_jobs=int(os.cpu_count()), backend="loky")(delayed(get_best_tblastn_result)(j, base_path, frst_evalue_thresold) for j in i)
    [selected_scaffolds.append(x) for x in temp]
del i, temp
selected_scaffolds = {k: v for d in selected_scaffolds for k, v in d.items()}
print("\nComplete\n")

#%% extracting scaffolds
# =============================================================================
# This for loop extracts scaffolds for ORFFinder
# =============================================================================
#% selected_scaffolds_new_seqs = pd.DataFrame(index=scaffolds.index, columns=scaffolds.columns)
print("\nExtracting scaffolds for ORFFinder\n")
selected_scaffolds_seqs  = {}
window_size = 5000
for i in tqdm(selected_scaffolds, unit = " sequences", position=0, leave=True):
    # print(f"Extracting a scaffold of window size {2*window_size} for hits for {i}")
    temp = pd.DataFrame(index=selected_scaffolds[i].index, columns = ["scaffold", "hit"])    
    
    for j in selected_scaffolds[i].index:        
        if selected_scaffolds[i]["subject id"][j] == "No hit found":
            temp["scaffold"][j] = "No hit found"
            temp["hit"][j] = "No hit found"
            continue
        key_1 = j.replace(f"{i}_", "")
        key_1 = key_1[:-4]
        seq = [x for x in genomes[f"{key_1}.fna"] if x.id in selected_scaffolds[i]["subject id"][j]]
        # input("j")  
        if type(selected_scaffolds[i]["subject id"][j]) is list:
            # input()
            scaffolds_for_selection = list(zip(selected_scaffolds[i]["subject id"][j], selected_scaffolds[i]["s. start"][j], selected_scaffolds[i]["s. end"][j]))
            temp1 = []
            for k in scaffolds_for_selection:
                
                if int(k[1]) < int(k[2]):
                    start = int(k[1])
                    end = int(k[2])
                if int(k[2]) < int(k[1]):
                    start = int(k[2])
                    end = int(k[1])
                    
                
                seq_n = [x for x in seq if x.id == k[0]]
                seq1 = str(seq_n[0].seq)
                
                seq_hit = seq1[start-1:end]
                    # seq1 = str(seq_n[0].seq)
              
                if len(seq1) < 2*window_size:
                    seq1 = seq1
                    temp1.append([seq1, seq_hit])
                    continue
                if start < window_size:
                    seq1 = seq1[0:start+window_size]
                    temp1.append([seq1, seq_hit])
                    continue
                if start + window_size > len(seq1):
                    seq1 = seq1[start-window_size:len(seq1)-1]
                    temp1.append([seq1, seq_hit])
                    continue
                
                seq1 = seq1[start-window_size:start+window_size]
                
                temp1.append([seq1, seq_hit])
                        
            temp["scaffold"][j] = [x[0] for x in temp1]
            temp["hit"][j] = [x[1] for x in temp1]
            continue
        
        
        if int(selected_scaffolds[i]["s. start"][j]) > int(selected_scaffolds[i]["s. end"][j]):
            start = int(selected_scaffolds[i]["s. end"][j])
            end = int(selected_scaffolds[i]["s. start"][j])
            
        else:
            start = int(selected_scaffolds[i]["s. start"][j])
            end = int(selected_scaffolds[i]["s. end"][j])
        
        seq1 = str(seq[0].seq)
        
        seq_hit = seq1[start-1:end]
        
        if len(seq1) < 2*window_size:
            # seq1 = seq1
            temp["scaffold"][j] = seq1
            temp["hit"][j] = seq_hit
            continue
        if start < window_size:
            seq1 = seq1[0:start+window_size]
            temp["scaffold"][j] = seq1
            temp["hit"][j] = seq_hit
            continue
        if start + window_size > len(seq1):
            seq1 = seq1[start-window_size:len(seq1)-1]
            temp["scaffold"][j] = seq1
            temp["hit"][j] = seq_hit
            continue
            
        
        seq1 = seq1[start-window_size:start+window_size]
        
        temp["scaffold"][j] = seq1
        temp["hit"][j] = seq_hit
        
    selected_scaffolds_seqs[i] = temp
del i, temp, j, key_1, seq, seq_hit, seq1, start, end, temp1, scaffolds_for_selection, k, seq_n, genomes
print("\nComplete\n")
#%% running ORFFinder
# =============================================================================
# This block of code uses ORFFinder to predict ORFs. The user must keep in mind
# that the ORF prediction itself will depend on the minimum length specified.
# To avoid setting the same ORF size for all proteins, minimum ORF length was
# specified based on the length of the tblastn hit. But, this seems to fail in
# several case (fail == ORF prediction gives 0 ORFs in a segment that is at
# least 6000 positions). For now, a standard minimum length of 75 is specified. 
# Even in this case, there are cases where no ORF is predicted. To counter 
# this, window size is increased from 3000 to 5000
# =============================================================================
scaffold_orfs = {}
print("\nStarting ORFFinder run\n")

with open(os.path.join(os.path.dirname(os.getcwd()), "orffinder_memory_usage.txt"), 'w') as f:
    f.write("Tracking memory usage for forward searches\n\n")
    f.write(f"{datetime.now()} - Current,Peak\n")
    f.close()

it = iter(selected_scaffolds_seqs)
new_dict = []
for i in range(0, len(selected_scaffolds_seqs), os.cpu_count()):
    new_dict.append({k:selected_scaffolds_seqs[k] for k in islice(it, os.cpu_count())})
del i, it

def orffinder_run(i, selected_scaffolds_seqs):
    temp = {}
    temp_orfs = {}
    # print(f"Running ORFFinder for {i}")
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        for j in selected_scaffolds_seqs[i].index:
            if selected_scaffolds_seqs[i]["scaffold"][j] == "No hit found":
                
                    temp[j] = "No hit found"
                    continue
            if type(selected_scaffolds_seqs[i]["scaffold"][j]) is list:
                temp1 = []
                
                for k in range(len(selected_scaffolds_seqs[i]["scaffold"][j])):                
                    seq = SeqRecord(Seq(selected_scaffolds_seqs[i]["scaffold"][j][k]))
                    temp1.append((orffinder.getORFProteins(seq, minimum_length=210, remove_nested=True, return_loci=True), (selected_scaffolds[i]["subject id"][j][k], selected_scaffolds[i]["s. start"][j][k], selected_scaffolds[i]["s. end"][j][k], len(selected_scaffolds_seqs[i]["scaffold"][j][k]))))
                    temp1[0]
                temp[j] = temp1
                continue
            seq = SeqRecord(Seq(selected_scaffolds_seqs[i]["scaffold"][j]))
            
            temp[j] = (orffinder.getORFProteins(seq, minimum_length=210, remove_nested=True, return_loci=True), (selected_scaffolds[i]["subject id"][j], selected_scaffolds[i]["s. start"][j], selected_scaffolds[i]["s. end"][j], len(selected_scaffolds_seqs[i]["scaffold"][j])))
    # print(f"Running ORFFinder for {i} - Complete")    
    temp_orfs[i] = temp
    return temp_orfs

tracemalloc.start()
scaffold_orfs = []
for i in tqdm(new_dict, unit = " sets", position=0, leave=True):    
    print(len(i))
    temp = Parallel(n_jobs=int(os.cpu_count()), backend='loky')(delayed(orffinder_run)(j, i) for j in i)
    [scaffold_orfs.append(x) for x in temp]
    with open(os.path.join(os.path.dirname(os.getcwd()), "orffinder_memory_usage.txt"), 'a') as f:
        f.write(f"\n{datetime.now()} - {round(tracemalloc.get_traced_memory()[0]/1e3, 4)} KB,{round(tracemalloc.get_traced_memory()[1]/1e3, 4)} KB\n")
        f.close()
del i, temp, new_dict
tracemalloc.stop()
scaffold_orfs = {k: v for d in scaffold_orfs for k, v in d.items()}
print("\nComplete\n")


#%% selecting orfs 
# =============================================================================
# ORFs are selected if they contain the midpoint of the selected frame and is
# the largest ORF in the set of predicted ORFs. If there are multiple hits 
# containing the midpoint of the selected frame, they will all be selected
# =============================================================================
print("\nSelecting predicted ORFs based on location of the mid-point\n")
scaffold_orfs_selected = {}
for i in tqdm(scaffold_orfs, unit = " sequences", position=0, leave=True):
    temp1 = {}
    
    for j in scaffold_orfs[i]:
        # input()
        # if j == "CBF1_YEAST_Jamesozyma_rosinii.txt":
        #     input("i")
        if scaffold_orfs[i][j] == "No hit found":
            temp1[j] = "No hit found"
            continue
        temp = []
        
        if type(scaffold_orfs[i][j]) is list:
            
            for k in scaffold_orfs[i][j]:            
                # if type(k) is list:
                for l in k[0]:
                    start = l["start"]
                    end = l["end"]
                    if start <= k[1][3]/2 <= end or end <= k[1][3]/2 <= start:
                        temp.append((l, k[1]))
                    # if k[1][3] == 2*window_size:                        
                    else:
                        if start <= k[1][3]/2 - len(selected_scaffolds_seqs[i]["hit"][j]) <= end or end <= k[1][3]/2 - len(selected_scaffolds_seqs[i]["hit"][j]) <= start:
                            temp.append((l, k[1]))
                        
                            
        if type(scaffold_orfs[i][j]) is not list:
            for k in scaffold_orfs[i][j][0]:
                start = k["start"]
                end = k["end"]
                
                
                if start <= scaffold_orfs[i][j][1][3]/2 <= end or end <= scaffold_orfs[i][j][1][3]/2 <= start:
                    temp.append((k, scaffold_orfs[i][j][1]))
                # if scaffold_orfs[i][j][1][3] == 2*window_size:
                else:
                    if start <= scaffold_orfs[i][j][1][3]/2 - len(selected_scaffolds_seqs[i]["hit"][j]) <= end or end <= scaffold_orfs[i][j][1][3]/2 - len(selected_scaffolds_seqs[i]["hit"][j]) <= start:
                        temp.append((k, scaffold_orfs[i][j][1]))
        if temp:
            temp1[j] = temp
        else:
            temp1[j] = "No hit found"
    scaffold_orfs_selected[i] = temp1
del i, j, k, temp, temp1, start, end, l
print("\nComplete\n")
#%% writing orfs out
os.chdir(os.path.dirname(base_path))
os.mkdir(f"{folder_name}_selected_orfs")
os.chdir(f"{folder_name}_selected_orfs")
temp = os.getcwd()

for i in scaffold_orfs_selected:
    os.mkdir(i)
    os.chdir(i)
    print(f"Writing selected ORFs for {i}\n")
    for j in scaffold_orfs_selected[i]:
        # input()
        if scaffold_orfs_selected[i][j] == "No hit found":
            # temp[j] = "No hit found"
            # os.chdir(temp)
            continue
        if not scaffold_orfs_selected[i][j]:
            # temp[j] = "No hit found"
            continue
        if len(scaffold_orfs_selected[i][j]) > 1:
            with open(f"{j}_orfs.fasta", 'w') as f:                
                p = 0
                completed = []
                for k in range(len(scaffold_orfs_selected[i][j])):                    
                    l = scaffold_orfs_selected[i][j][k]
                    if l[0]["protein"] in completed:
                        continue
                    else:
                        # input()
                        SeqIO.write(SeqRecord(l[0]["protein"], id=f"{j}_{p}", description=f"[From tblastn search results - scaffold={l[1][0]}, start={l[1][1]}, end={l[1][2]}, orf_start={l[0]['start']}, orf_end={l[0]['end']}, (predicted orf centered around start. start coordinate is equivalent of position {window_size})]"), f, 'fasta')
                        completed.append(l[0]["protein"])
                        p += 1
                f.close()
                # os.chdir(temp)
                continue
        with open(f"{j}_orfs.fasta", 'w') as f:
            p = 0
            completed = []
            for k in range(len(scaffold_orfs_selected[i][j])):
                if scaffold_orfs_selected[i][j][k][0]["protein"] in completed:
                    continue
                else:
                    # input()
                    SeqIO.write(SeqRecord(scaffold_orfs_selected[i][j][k][0]["protein"], id=f"{j}_{p}", description=f"[From tblastn search results - scaffold={scaffold_orfs_selected[i][j][k][1][0]}, start={scaffold_orfs_selected[i][j][k][1][1]}, end={scaffold_orfs_selected[i][j][k][1][2]}, orf_start={scaffold_orfs_selected[i][j][k][0]['start']}, orf_end={scaffold_orfs_selected[i][j][k][0]['end']} (predicted orf centered around start. start coordinate is equivalent of position {window_size})]"), f, 'fasta')
                    completed.append(scaffold_orfs_selected[i][j][k][0]["protein"])
                    p += 1
            f.close()
    os.chdir(temp)
#%% writing out all orfs
os.chdir(os.path.dirname(base_path))  
os.mkdir(f"{folder_name}_all_orfs")
os.chdir(f"{folder_name}_all_orfs")
temp = os.getcwd()
#%
for i in scaffold_orfs:
    os.mkdir(i)
    os.chdir(i)
    print(f"Writing all ORFs for {i}\n")
    for j in scaffold_orfs[i]:
        if scaffold_orfs[i][j] == "No hit found":
            # temp[j] = "No hit found"
            # os.chdir(temp)
            continue
        if not scaffold_orfs[i][j]:
            # temp[j] = "No hit found"
            continue
        if type(scaffold_orfs[i][j]) is list:
            with open(f"{j}_orfs.fasta", 'w') as f:                
                p = 0
                completed = []
                for k in range(len(scaffold_orfs[i][j])):                    
                    for l in scaffold_orfs[i][j][k][0]:
                        if l["protein"] in completed:
                            continue
                        else:
                            # input()
                            SeqIO.write(SeqRecord(l["protein"], id=f"{j}_{p}", description=f"[From tblastn search results - scaffold={scaffold_orfs[i][j][k][1][0]}, start={scaffold_orfs[i][j][k][1][1]}, end={scaffold_orfs[i][j][k][1][2]}, orf_start={l['start']}, orf_end={l['end']} (predicted orf from scaffold sequence of length {scaffold_orfs[i][j][k][1][3]})]"), f, 'fasta')
                            completed.append(l["protein"])
                            p += 1
                f.close()
                # os.chdir(temp)
                continue
        with open(f"{j}_orfs.fasta", 'w') as f:
            p = 0
            completed = []
            for k in scaffold_orfs[i][j][0]:
                if k["protein"] in completed:
                    continue
                else:
                    # input()
                    SeqIO.write(SeqRecord(k["protein"], id=f"{j}_{p}", description=f"[From tblastn search results - scaffold={scaffold_orfs[i][j][1][0]}, start={scaffold_orfs[i][j][1][1]}, end={scaffold_orfs[i][j][1][2]}, orf_start={k['start']}, orf_end={k['end']} (predicted orf from scaffold sequence of length {scaffold_orfs[i][j][1][3]})]"), f, 'fasta')
                    completed.append(k["protein"])
                    p += 1
            f.close()
    os.chdir(temp)
del i, j, k, f, p, l, temp, completed
#%% writing a log file

os.chdir(os.path.dirname(base_path)) 
with open(rf"{folder_name}_log_{datetime.today().strftime('%d_%m_%Y')}.txt", 'w') as f:
    f.write("*****Log file*****\n\n")
    f.write("Analysis of tblastn search results\n\n")
    f.write(f"tblastn searches were used to identify the scaffold(s) in the genome assembl(y/ies) containing homolog(s) for the gene(s) of interest. Within the scaffold, a segment of length {2*window_size} positions centered on the starting cooardinate of the best-scoring hits is extracted. ORFFinder (a python version of the NCBI tool) was used to identify ORFs within the selected segment. Only ORFs containing the midpoint of the segment were selected and written out for further analysis (reciprocal BLAST).\n\n")
    f.write(f"****Species used - {len(list_of_genomes.keys())}****\n\n")
    [f.write(f"{x[:-4].replace('_', ' ')}\n") for x in list(list_of_genomes.keys())]
    f.write(f"\n\n****The {len(base_list)} genes of interest in this analysis were selected from Saccharomyces cerevisiae****\n\n")
    [f.write(f'{x}\n') for x in base_list]
    f.write("\n\n****The selected ORFs were written out****\n\n")
    f.close()
