# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 12:59:54 2024

@author: Kausthubh R
"""

import os
import subprocess as sp
from concurrent.futures import ThreadPoolExecutor as tpe
from datetime import datetime
import tracemalloc
from tqdm import tqdm

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
        
    return futures

# =============================================================================
# It is recommended to commence your run in a fresh directory which contains 
# only one subdirectory containing fasta files (one sequence per fasta file)
# for your proteins of interest
# =============================================================================

path2markerseqs = input("Please enter the path to the directory containing your marker sequences:\n")
list_of_files = [os.path.join(path2markerseqs, x) for x in os.listdir(path2markerseqs)]

os.chdir(os.path.dirname(os.getcwd()))
os.mkdir("phmmer_searches")
os.chdir("phmmer_searches")

path2db = input("Please enter the path to the BLAST database against which you want to perform phmmer searches:\n")

if os.name == 'nt':
    cmd_list = [f'wsl phmmer --tblout {x[x.rfind("markers")+8:]}.txt {x} {path2db}' for x in list_of_files]
else:
    cmd_list = [f'phmmer --tblout {x[x.rfind("/")+1:]}.txt {x} {path2db}' for x in list_of_files]
    
    
with open(os.path.join(os.path.dirname(os.getcwd()), "phmmer_searches_memory_usage.txt"), 'w') as f:
    f.write("Tracking memory usage for forward searches\n\n")
    f.write(f"{datetime.now()} - Current,Peak\n")
    f.close()
    #%%
# =============================================================================
# tracemalloc is used to monitor memory usage to enable easier debugging if
# the parallelized phmmer searches fail.
# =============================================================================
tracemalloc.start()
no_of_workers = os.cpu_count()
# Parallel(n_jobs=16)(delayed(sp.getoutput)(i) for i in cmd_list)
cmd_list = [cmd_list[i:i + 100] for i in range(0, len(cmd_list), 100)]
for i in tqdm(cmd_list, unit=" commands"):
    pull_run(no_of_workers, i)    
    with open(os.path.join(os.path.dirname(os.getcwd()), "phmmer_searches_memory_usage.txt"), 'a') as f:
        # old = f.read()
        f.write(f"\n{datetime.now()} - {round(tracemalloc.get_traced_memory()[0]/1e3, 4)} KB,{round(tracemalloc.get_traced_memory()[1]/1e3, 4)} KB\n")
        f.close()
del i, f

tracemalloc.stop()


