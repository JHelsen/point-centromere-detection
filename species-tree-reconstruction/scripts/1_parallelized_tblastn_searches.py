# -*- coding: utf-8 -*-
"""

@author: Kausthubh R
"""
# =============================================================================
# This script performs tblastn searches for a list of protein sequences against 
# a list of species' genomes.
# 
# The species need to be organized as a list of BLAST databases in a separate
# folder "genomes_blast_dbs", one database per species.
# 
# The sequences need to be in a separate folder "genes_of_interest_sequences" 
# which should contain nothing but the sequence files as FASTA files, one file 
# per sequence.
# =============================================================================
import subprocess as sp
import os
from concurrent.futures import ThreadPoolExecutor as tpe

#% function cell
def exec_(cmd:str):
    '''
    single command run
    command must be prepared as subprocess.call
    '''
    sp.getoutput(cmd)

        
def pull_run(obj1:int, obj2:list):
    print(f"Number of searches being performed: {len(obj2)}\n")
    if len(obj2) > obj1:
        with tpe(max_workers=obj1) as exe:
            futures = exe.map(exec_, obj2)
            # wait(futures)
    else:
        with tpe(max_workers=len(obj2)) as exe:
            futures = exe.map(exec_, obj2)
        
    return futures

# =============================================================================
# Getting paths to the BLAST databases
# =============================================================================
start_point = os.getcwd()

path2blastdbs = input("Please enter the path to the directory containing your BLAST databases:\n")

unique_dbs = list(dict.fromkeys([x[:x.rfind(".")] for x in os.listdir(path2blastdbs)]))
list_of_blast_dbs = {x : os.path.join(path2blastdbs, x) for x in unique_dbs}
del unique_dbs

list_of_blast_dbs = os.listdir() 


os.chdir(start_point)


# =============================================================================
# It is recommended to commence your run in a fresh directory which contains 
# only one subdirectory containing fasta files (one sequence per fasta file)
# for your proteins of interest
# =============================================================================
# os.chdir(r"orthogroups_12092024")
path2seqs = input("Please enter the path to the directory containing your query sequences:\n")
run_title = input("Please provide a relevant tile for your tblastn run:\n")
os.chdir(path2seqs)
#

marker_seqs = [os.path.join(os.getcwd(), x) for x in os.listdir()]

os.chdir(os.path.dirname(os.getcwd()))

#%
os.mkdir(rf"{run_title}_tblastn_results")
os.chdir(rf"{run_title}_tblastn_results")


base_path = os.getcwd()

# =============================================================================
# Performing parallelized tblastn searches using sp.getoutput()
# =============================================================================
for i in marker_seqs:
    name = i[i.rfind("/")+1:]
    name = name[:name.find(".")]
    print(name)
    os.mkdir(name)
    os.chdir(name)
    cmd_list = [f"tblastn -db {list_of_blast_dbs[x]} -query {i} -outfmt '7 evalue bitscore score length pident nident qstart qend sstart send sseqid sseq' -out {name}_{x[:x.rfind('.')]}.txt" for x in list_of_blast_dbs]
    pull_run(os.cpu_count(), cmd_list)
    os.chdir(base_path)



