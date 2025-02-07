# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 21:41:17 2024

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

# =============================================================================
# It is recommended to commence your run in the same directory which contains 
# the results of the 1_parallelized_tblastn_searches.py script and the 
# 2_tblastn_results_analysis_orf_extraction.py script
# =============================================================================

cerevisiae_blast_db = input("Please enter the path to the Saccharomyces cerevisiae proteome used for reciprocal BLASTp:\n")
path2selectedorfs = input("Please enter the path to the directory containing ORFs extracted using ORFFinder:\n")
os.chdir(path2selectedorfs)
folder_name = os.path.basename(os.getcwd())

#%
marker_seqs = {}

for i in [x for x in os.listdir() if not x.startswith(".")]:
    marker_seqs[i] = [os.path.join(os.path.join(os.getcwd(), i), x) for x in os.listdir(os.path.join(os.getcwd(), i)) if not x.startswith(".")]
   
del i
#%

os.chdir(os.path.dirname(os.getcwd()))

os.mkdir(rf"{folder_name}_reciprocal_blast")
os.chdir(rf"{folder_name}_reciprocal_blast")

base_path = os.getcwd()

for i in marker_seqs:
    
    print(f"Performing reciprocal blast for {i}")
    os.mkdir(i)
    os.chdir(i)
 
    cmd_list = [f"blastp -db {cerevisiae_blast_db} -query {x} -outfmt '7 evalue bitscore score length pident nident qstart qend sstart send sseqid sseq' -out {x[x.rfind('/')+1:-15]}_reciprocal_blast.txt" for x in marker_seqs[i]] # for marker genes used for species tree reconstruction
    
    pull_run(os.cpu_count(), cmd_list)
    print(f"Completed reciprocal blast for {i}")
    
    os.chdir(base_path)

