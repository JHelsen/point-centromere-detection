# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 13:29:16 2024

@author: Kausthubh R
"""

#%%
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
from orffinder import orffinder
from tqdm import tqdm
from joblib import Parallel, delayed
import warnings
warnings.filterwarnings("ignore")


# =============================================================================
# This block of code reads tblastn search results [command used - tblastn -db 
# {list_of_blast_dbs[j]} -query {i} -outfmt '7 evalue bitscore score length 
# pident nident qstart qend sstart send sseqid sseq' -out {name}_{j}.txt] and 
# makes a dictionary of type Protein:Dataframe
# =============================================================================

path2genomes = input("Please enter the path to the directory containing your genomes for the analysis:\n")
list_of_genomes = dict([(x, os.path.join(path2genomes, x)) for x in os.listdir(path2genomes)])
print("\nLoading genomes\n")
genomes = {}
for i in tqdm(list_of_genomes, unit=" genomes", position=0, leave=True):
    genomes[i] = [x for x in SeqIO.parse(list_of_genomes[i], 'fasta')]
    # print(i)
del i

# genomes = list(genomes.keys())
print("\nComplete\n")
#%%
# =============================================================================
# Importing protein sequences
# =============================================================================
path2proteinseqs = input("Please enter the path to the directory containing the protein homologs from the ORFFinder extraction after selecting the correct, reciprocally mapped ones:\n")
os.chdir(path2proteinseqs)
#%
curated_seqs = {}
for i in tqdm([x for x in os.listdir() if not x.startswith(".")]):
    curated_seqs[i] = [{x[:x.rfind(".")] : y for y in SeqIO.parse(os.path.join(os.path.join(os.getcwd(), i), x), 'fasta')} for x in os.listdir(os.path.join(os.getcwd(), i)) if not x.startswith(".")]
del i
curated_seqs = {x : {k: v for d in curated_seqs[x] for k, v in d.items()} for x in curated_seqs}

os.chdir(os.path.dirname(os.getcwd()))
#%% importing tblastn results
# =============================================================================
# This block of code only travels to each directory, opens the tblastn results
# in a dataframe, and only retains the hit with the lowest evalue 
# =============================================================================

os.chdir(os.path.dirname(os.getcwd()))
run_title = os.path.basename(os.getcwd())

print("\nImporting tblastn results\n")
os.chdir(f"{run_title}_seqs_reciprocal_blast") # ONLY for manually curated genes of interest
selected_scaffolds = {}
base_list = [x for x in os.listdir() if not x.startswith(".")]
base_list_import = [base_list[i:i + int(os.cpu_count())] for i in range(0, len(base_list), int(os.cpu_count()))]
base_path = os.getcwd() #the starting folder should be the folder containing the tblastn search results
def get_best_tblastn_result(i, base_path):
    os.chdir(i)
    temp_dict= {}
    # print(f"Importing scaffolds for {i}\n")
    temp2 = pd.DataFrame(columns=["subject id", "s. start", "s. end", "evalue"], index = os.listdir())
    # temp2 = []
    for j in os.listdir():
        # input()
        file = open(j, 'r').read().split("\n")
        if len(file) <= 6:
            temp2["subject id"][j] = "No hit found"            
            continue
        
        headers = file[3]
        headers = headers[headers.find(":")+1:].split(", ")
        headers[0] = "evalue"
        
        values = file[5:len(file)-2]
        values = [x.split("\t") for x in values]
        values = [list(map(float, [x[y] for y in range(len(x)) if not y > 9])) + [x[y] for y in range(len(x)) if y > 9]  for x in values]
        results = pd.DataFrame(values, columns = headers)
        temp_2 = results[results["evalue"] == min(list(results["evalue"]))]
        # input()
        # temp = results[results["evalue"] < frst_evalue_thresold]
        temp2["subject id"][j] = temp_2["subject id"][temp_2.index[0]]
        temp2["s. start"][j] = int(temp_2["s. start"][temp_2.index[0]])
        temp2["s. end"][j] = int(temp_2["s. end"][temp_2.index[0]])
        temp2["evalue"][j] = float(temp_2["evalue"][temp_2.index[0]])
            
    temp_dict[i] = temp2
    os.chdir(base_path)
    return temp_dict

selected_scaffolds = []
# p = 0
for i in tqdm(base_list_import, unit=" sequences", position=0, leave=True):
    temp = Parallel(n_jobs=int(os.cpu_count()), backend="loky")(delayed(get_best_tblastn_result)(j, base_path) for j in i)
    [selected_scaffolds.append(x) for x in temp]
    # print(p)
    # p += 1
del i, temp
selected_scaffolds = {k: v for d in selected_scaffolds for k, v in d.items()}
print("\nComplete\n")
os.chdir(os.path.dirname(os.getcwd()))
#%% extracting scaffolds
# =============================================================================
# This for loop extracts scaffolds for ORFFinder. It first checks if the 
# tblastn hit matches the entire protein when translated. If it does not, it
# then runs ORFFinder
# =============================================================================
#% selected_scaffolds_new_seqs = pd.DataFrame(index=scaffolds.index, columns=scaffolds.columns)
print("\nRunning ORFFinder for cases where the hit doesn't match the entire query sequence\n")
selected_orf_seqs  = {}
window_size = 3000
for i in tqdm(selected_scaffolds, unit = " sequences", position=0, leave=True):
    # input()
    # print(f"Extracting a scaffold of window size {2*window_size} for hits for {i}")
    temp = pd.DataFrame(index=selected_scaffolds[i].index, columns = ["orf", "reversed"])    
    
    for j in selected_scaffolds[i].index:  
        # input()
        # if j == "CBF1_AA_Zygosaccharomyces_parabailii_2.txt":
        #     input()
        if selected_scaffolds[i]["subject id"][j] == "No hit found":
            temp["orf"][j] = "No hit found"
            temp["reversed"][j] = "No hit found"
            continue
        key_1 = j.replace(f"{i}_", "")
        og = key_1[:key_1.rfind(".")]
        key_1 = key_1[:key_1.rfind(".")]
        # key_1 = key_1.replace("Fragment", "")
        
        if key_1[-2] == "_":
            key_1 = key_1[:key_1.rfind("_")]
            
        
        seq = [x for x in genomes[f"{key_1}.fna"] if x.id in selected_scaffolds[i]["subject id"][j]]
        # input("j")  
        
        reverse_flag = False
        if int(selected_scaffolds[i]["s. start"][j]) > int(selected_scaffolds[i]["s. end"][j]):
            start = int(selected_scaffolds[i]["s. end"][j])
            end = int(selected_scaffolds[i]["s. start"][j])
            reverse_flag = True
        else:
            start = int(selected_scaffolds[i]["s. start"][j])
            end = int(selected_scaffolds[i]["s. end"][j])
        # input()
        
        seq1 = str(seq[0].seq)
        
        seq_hit = seq1[start-1:end]

        # =============================================================================
        # Checking if the extracted ORF matches the query protein
        # =============================================================================        
        query_protein = str(curated_seqs[i][[x for x in curated_seqs[i] if x.startswith(og)][0]].seq)
        
        if reverse_flag:
            seq_hit_protein = str(Seq(seq_hit).reverse_complement().translate())
        else:
            seq_hit_protein = str(Seq(seq_hit).translate())
        # input()
        if query_protein == seq_hit_protein:
            if reverse_flag:
                temp["orf"][j] = SeqRecord(Seq(seq_hit).reverse_complement())
                temp["reversed"][j] = "Y"
            else:
                temp["orf"][j] = SeqRecord(Seq(seq_hit))
                temp["reversed"][j] = "N"
        else:
            if len(seq1) < 2*window_size:
                scaffold = seq1
            
            elif start < window_size:
                scaffold = seq1[0:start+window_size]
               
            elif start + window_size > len(seq1):
                scaffold = seq1[start-window_size:len(seq1)-1]
            
            else:
                scaffold = seq1[start-window_size:start+window_size]
                
            scaffold = SeqRecord(Seq(scaffold))
            temp_orfs = orffinder.getORFProteins(scaffold, minimum_length=len(query_protein), remove_nested=True, return_loci=True)
            
            orffinder_change = False
            for k in temp_orfs:
                for_comparison = str(k["protein"])
                for_comparison = "".join([x for x in for_comparison if x.isalpha()])                
                # input()
                if for_comparison == query_protein:                    
                    if reverse_flag:
                        seq_hit = str(scaffold.seq)[k["end"]:k["start"]-1]
                        orffinder_change = True
                    else:
                        seq_hit = str(scaffold.seq)[k["start"]-1:k["end"]]
                        orffinder_change = True
            
            if not orffinder_change:
                temp["orf"][j] = "ORFFinder failed to recover any ORF that perfectly matches the query protein"
                temp["reversed"][j] = "N"
            
            if orffinder_change:
                if reverse_flag:
                    if "".join([x for x in str(Seq(seq_hit).reverse_complement().translate()) if x.isalpha()]) == query_protein:
                        temp["orf"][j] = SeqRecord(Seq(seq_hit).reverse_complement())
                        temp["reversed"][j] = "Y"
                    else:
                        temp["orf"][j] = "ORFFinder has not recovered any ORF that perfectly matches the query protein"
                        temp["reversed"][j] = "Y"
                else:
                    if "".join([x for x in str(Seq(seq_hit).translate()) if x.isalpha()]) == query_protein:
                        temp["orf"][j] = SeqRecord(Seq(seq_hit))
                        temp["reversed"][j] = "N"
                    else:
                        temp["orf"][j] = "ORFFinder has not recovered any ORF that perfectly matches the query protein"
                        temp["reversed"][j] = "N"
    selected_orf_seqs[i] = temp
del i, temp, j, key_1, seq, seq_hit, seq1, start, end, k, 
print("\nComplete\n")
#%%
with open(r"ORF_recovery_log.txt", 'w') as f:
    f.write("ORFs could not be identified for the following sequences in each gene of interest\n\n")
    f.close()
    
path2log = os.path.join(os.getcwd(), r"ORF_recovery_log.txt")

os.mkdir(f"{run_title}_seqs_orfs")
os.chdir(f"{run_title}_seqs_orfs")

base_path = os.getcwd()

for i in selected_orf_seqs:
    os.mkdir(i)
    os.chdir(i)
    with open(path2log, 'a') as g:
        g.write(f"{i}\n")
        g.close()
    with open(f"{i}_orfs.fasta", 'w') as f:
        for j in selected_orf_seqs[i].index:
            # input()
            if type(selected_orf_seqs[i]["orf"][j]) is str:
                with open(path2log, 'a') as g:
                    g.write(f"{j}\n")
                    g.close()
                    continue
            selected_orf_seqs[i]["orf"][j].id = j[:j.rfind(".")]
            if selected_orf_seqs[i]["reversed"][j] == "Y":
                selected_orf_seqs[i]["orf"][j].description = "orf from reverse strand"
            else:
                selected_orf_seqs[i]["orf"][j].description = ""
            try:
                SeqIO.write(selected_orf_seqs[i]["orf"][j], f, 'fasta-2line')
            except TypeError:
                with open(path2log, 'a') as g:
                    g.write(f"{j}\n")
                    g.close()
    with open(path2log, 'a') as g:
        g.write("\n\n")
        g.close()
    os.chdir(base_path)      
                

