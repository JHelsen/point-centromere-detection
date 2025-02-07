# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 01:20:12 2024

@author: Kausthubh R
"""

import subprocess as sp
import os
from concurrent.futures import ThreadPoolExecutor as tpe


#% function cell
def exec_(cmd):
    '''
    single command run
    command must be prepared as subprocess.call
    '''
    sp.getoutput(cmd)

        
def pull_run(obj1, obj2):
    print(f"Number of searches being performed: {len(obj2)}\n")
    if len(obj2) > obj1:
        with tpe(max_workers=obj1) as exe:
            futures = exe.map(exec_, obj2)
            # wait(futures)
    else:
        with tpe(max_workers=len(obj2)) as exe:
            futures = exe.map(exec_, obj2)
    
    del futures



start_point = os.getcwd()

# blast_dbs = [os.path.join(r"Z:\Kausthubh\1_postdocs\Jana\Centromere_project_phylogeny\Genomes", x) for x in os.listdir(r"Z:\Kausthubh\1_postdocs\Jana\Centromere_project_phylogeny\Genomes")]
blast_dbs = [os.path.join(r"/g/dey/Kausthubh/1_postdocs/Jana/Centromere_project_phylogeny/blast_dbs", x) for x in os.listdir(r"/g/dey/Kausthubh/1_postdocs/Jana/Centromere_project_phylogeny/Genomes")]
folder_name = r"curated_seqs"
os.chdir(rf"{folder_name}_split")
#%%
curated_seqs = {}
for i in [x for x in os.listdir() if not x.startswith(".")]:
    curated_seqs[i] = [os.path.join(os.path.join(os.getcwd(), i), x) for x in os.listdir(os.path.join(os.getcwd(), i)) if not x.startswith(".")]
del i

os.chdir(start_point)
#%
os.mkdir(rf"{folder_name}_seqs_reciprocal_blast")
os.chdir(rf"{folder_name}_seqs_reciprocal_blast")

base_path = os.getcwd()
# print(base_path)
for i in curated_seqs:
    print(i)
    species_name = [x[x.rfind("/")+1:] for x in curated_seqs[i]]
    species_name = [x[:x.rfind(".")] for x in species_name]
    species_name = [x[:x.rfind("_")] if "Fragment" in x else x for x in species_name]
    species_name = [x[:x.rfind("_")] if "fragment" in x else x for x in species_name]
    species_name = [x[:x.rfind("_")] if "fromOtherAssembly" in x else x for x in species_name]
    species_name = [x[:x.rfind("_")] if "differentAssembly" in x else x for x in species_name]
    species_name = [x[:x.rfind("_")] if "OtherAssembly" in x else x for x in species_name]
    species_name = [x[:x.rfind("_")] if "PossibleFalsePos" in x else x for x in species_name]
    species_name = [x[:x.rfind("_")] if "Broken" in x else x for x in species_name]
    # species_list = [x if x.count("_") == 1 else x[:x.rfind("_")] for x in species_name]
    # species_list = list(set(species_list))
    
    
    # input()
    print(f"Performing reciprocal blast for {i}")
    try:
        os.mkdir(i)
    except FileExistsError:
        print(f"tblastn_searches_completed for {i}")
        os.chdir(base_path)
        continue
    
    os.chdir(i)
    cmd_list = []
    for j in species_name:
        temp = [x for x in curated_seqs[i] if f"{j}." in x or j in x][0]
        temp_name = [x for x in species_name if x in temp][0]
        # if j == "Torulaspora_sp2017_383":
        #     input()
        while temp_name[-2] == "_":
            temp_name = temp_name[:temp_name.rfind("_")]
       
        cmd_list.append(f"tblastn -db {[x for x in blast_dbs if temp_name in x][0]} -query {temp} -outfmt '7 evalue bitscore score length pident nident qstart qend sstart send sseqid sseq' -out {i}_{j}.txt")
    pull_run(os.cpu_count(), cmd_list)
    print(f"Completed reciprocal blast for {i}\n")
    
    os.chdir(base_path)
