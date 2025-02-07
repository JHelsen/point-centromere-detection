# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 22:06:12 2024

@author: Kausthubh R
"""


import os
import subprocess as sp
from Bio import SeqIO

trimv = sp.getoutput("trimal --version")
folder_path = os.getcwd()

align_check = [x for x in os.listdir() if x.endswith("_trim.phy")]
if align_check:
    with open("error.txt", 'w') as f:
        f.write("this folder already seems to contain alignments. Please check")
        f.close()
    exit()
    
name_list = [f for f in os.listdir() if f.endswith('.fasta')]

if not name_list:
    with open("error.txt", 'w') as f:
        f.write("this folder does not seem to have any fasta files. Please check")
        f.close()
    exit()

name = name_list[0]

if not name:
    exit()
name_out = name[:-6] + "_align.phy"
print(name)
temp_record = []
for record in SeqIO.parse(open(name), 'fasta'):
    temp_record.append(record)
del record

# =============================================================================
# Exiting if the number of sequences available is less than 4
# =============================================================================

if len(temp_record) <= 4:
    with open("{0}_record.txt".format(name), 'w') as f:
        f.write("\nThe number of available sequences is less than 4. It is therefore meaningless to build alignments / trees with bootstrap\n")
    exit()

temp_record = list(filter(None, temp_record))

temp_record_len = []
for i in temp_record:
    temp_record_len.append(len(i))
del i


# =============================================================================
# Removing sequences that are smaller than 10% of the average length of the 
# sequences in the dataset
# =============================================================================
avg = sum(temp_record_len)/(len(temp_record_len))
print(avg)
tiny_seqs = []
y = 0.1
for i in range(len(temp_record)):
    if len(temp_record[i]) < y*avg:
        tiny_seqs.append(temp_record[i])
        temp_record[i] = []
del i
temp_record = list(filter(None, temp_record))

if tiny_seqs:
    with open(f"{name[:name.rfind('.')]}_small_seqs.fasta", "w") as f:
        [SeqIO.write(x, f, 'fasta') for x in tiny_seqs]
        f.close()
        

if tiny_seqs:
    name_updated = f"{name[:name.rfind('.')]}_updated.fasta"

    with open((name_updated), 'w') as f:
        for i in temp_record:
            SeqIO.write(i, f, 'fasta')
        f. close()
    del f, i

else:
    name_updated = name


updated_len = len(temp_record)
print("Aligning the updated {0} sequences\n".format(name))

name_out = name_updated[:name_updated.rfind(".")] + "_align.phy"

# =============================================================================
# Exiting if the number of sequences available is more than 3000
# =============================================================================

if len(temp_record) > 3000:
    with open("{0}_record.txt".format(name), 'w') as f:
        f.write("\nThe number of available sequences is greater than 3000. It is therefore too computationally intensive to build alignments / trees with bootstrap\n")
    exit()


#mafft
if len(temp_record) <= 500:
    mafft_cmd = ("mafft --maxiterate 1000 --bl 45 --genafpair --reorder {0} > {1}").format(name_updated, name_out)
if len(temp_record) >= 500:
    mafft_cmd = ("mafft --retree 2 --maxiterate 1000 --bl 45 --reorder {0} > {1}").format(name_updated, name_out)

x = sp.getoutput(mafft_cmd) 

with open("mafft_record_{0}.txt".format(name_out[:-4]), 'w') as f:
    f.write(x)
    f.close()
del f, x

print("Trimming the genafpair aligned {0} sequences\n".format(name))
name_trim = name_out[:-4]+"_trim"

trim_cmd = "trimal -in {0} -out {1}.phy -gappyout -htmlout {1}_log.html".format(name_out, name_trim)
sp.getoutput(trim_cmd)

# =============================================================================
# Exiting if the number of sequences available is less than 4
# =============================================================================
if len(temp_record) < 4:
    with open("{0}_record.txt".format(name), 'w') as f:
        f.write("\nThe number of available sequences is less than 4. It is therefore meaningless to try and build trees with bootstrap\n")
    exit()

# =============================================================================
# Redoing the trimming if >70% of positions have been removed
# =============================================================================
retention_proportion = 0.3
retention_flag = False
print(f"\nChecking if trimmed MSA retains at least {retention_proportion*100}% of the positions from untrimmed MSA\n")
og_length = [len(x) for x in SeqIO.parse(name_out, 'fasta')][0]
new_length = [len(x) for x in SeqIO.parse("{0}.phy".format(name_trim), 'fasta')][0]
if new_length >= retention_proportion*og_length:
    print(f"Trimmed MSA retains at least {retention_proportion*100}% of the positions from untrimmed MSA\n")

if new_length < retention_proportion*og_length:
    retention_flag = True
    gap_threshold = 0.4
    print(f"Trimmed MSA does not retain at least {retention_proportion*100}% of the positions from untrimmed MSA\nRetrimming with a gap threshold parameter -gt ")
    os.remove("{0}.phy".format(name_trim))
    os.remove("{0}_log.html".format(name_trim))
    name_trim = name_out[:-4]+"_gt_trim"
    trim_cmd = "trimal -in {0} -out {1}.phy -htmlout {1}_log.html -gt {2}".format(name_out, name_trim, gap_threshold)
    sp.getoutput(trim_cmd)

# =============================================================================
# Tree building
# =============================================================================
print("Building trees using gappyout trimmed alignments of the {0} sequences\n\n\n".format(name))


tree_cmd = ('FastTree -spr 4 -mlacc 2 -slownni -n 1000 -gamma -log {0}_fasttree_log.txt {0}.phy > {0}_fasttree.tree'.format(name_trim))
sp.getoutput(tree_cmd)

with open("{0}_record.txt".format(name[:name.rfind('.')]), 'w') as f:
    f.write("Record of commands used for aligning, trimming, and building trees for {0}\n\n\n".format(name))
    
    f.write("Alignment tool used: MAFFT (version specified in tool's log file mafft_record_{0}.txt)\n".format(name_out[:-4]))
    f.write("Alignment tool command used: {0}\n".format(mafft_cmd))
    
    f.write("Aligment trimming tool used: TrimAl (version specified in tool's HTML output) \n")
    f.write("Aligment trimming tool used: {0}) \n".format(trimv))
    f.write("Alignment trimming tool command used: {0}\n".format(trim_cmd))
    
    if retention_flag:
        f.write("While using TrimAl, the gappyout method removed too many positions. So, a gap threshold (gt) of {0} was used\n".format(gap_threshold))
    
    f.write("Tree building tool used: FastTree (version specified in tool's log file)\n")
    f.write("Tree building tool command used: {0}\n".format(tree_cmd))
    f.write("Average protein length {0}\n\n\n".format(avg))
    f.write('Proteins that were smaller than {0} of the average protein length were removed from the alignment\n\n\n'.format(y))
    
    if tiny_seqs:
        f.write(f"The following sequences were removed prior to building alignments as the proteins were smaller than {y} times the average length ({avg} amino acids) of the proteins in the set: \n\n")
        [f.write(f"{x.id}\n") for x in tiny_seqs]
    
    f.close()
