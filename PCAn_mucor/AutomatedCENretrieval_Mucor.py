#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:33:28 2024

@author: janahelsen
"""

###########
# Libraries #
###########

import os
from Bio import SeqIO
import pandas as pd


##############################
# PARAMETERS AND INPUT FILES #
##############################

Genome = input("\nEnter input path to genome assembly file here. Please include the file name with extension in the path\n")

MotifsAndThresholds = {'Mucor motifs all':["UsedMotifs/250528_Mucormotif41bp_all.meme"],
                       'Mucor circinelloides (for 2)':["UsedMotifs/Mcircmotif41bp_2.meme"],
                       'Mucor circinelloides (other 3)':["UsedMotifs/Mlucmotif41bp_other3.meme"],
                       'Mucor racemosus':["UsedMotifs/Mracemosus41bp.meme"],
                       'Parasitella chaetocladium':["UsedMotifs/ParasitellaChaetocladium41bp.meme"],
                       }

available = list(MotifsAndThresholds.keys())
for i in range(len(available)):
    print(f"\n{i} - {available[i]}")
del i

Motif = MotifsAndThresholds[available[int(input("\nPlease select an option from the above:\n"))]][0]


# How many bp upstream to extract
Motifup = 0 #default 750

# How many bp downstream of CDEI to extract
Motifdn = 750 #default 750

#####################################################
# CANDIDATE CEN EXTRACTION: FIMO, DEFINING FUNCTION #
#####################################################

def FIMO2(Genome, Motif, Motifup, Motifdn):

    command1 = "fimo --thresh 1.0E-8 --oc FIMOresults "+ Motif + " "+ Genome
    os.system(command1)
    
    fimosearchResult = 'FIMOresults/fimo.txt'
    Queries = pd.read_csv(fimosearchResult, skiprows=0, delimiter="\t")
    
    Sequences = []
    Lines = []
    ListSeq = []
    
    for index, row in Queries.iterrows():
        ID = row[1]
        strand = row[4]
        start = row[2]
        end = row[3]
        QualityScoreI = row[5]
        if strand == '+':
            Newstart = start - Motifup - 1
            Newend = end + Motifdn
        if strand == '-':
            Newstart = start -1 - Motifdn
            Newend = end + Motifup
        if Newstart < 0:
            Newstart = 0
        for seq in SeqIO.parse(Genome,"fasta"):
            if seq.id == ID:
                if strand == '+':
                    sequence = str(seq.seq[Newstart:Newend])
                if strand == '-':
                    seqRC = seq.seq[Newstart:Newend]
                    sequence = str(seqRC.reverse_complement())          
                Sequences.append(sequence)
        Lines.append(">"+ ID + "_" + str(index+1) + "\n" + sequence)
        ListSeq.append([ID + "_" + str(index+1), Newstart, Newend, sequence, QualityScoreI])
    
    with open('SequenceHits.fsa','a') as outfile:
        outfile.writelines(i + '\n' for i in Lines)
        outfile.close
       
    ColumnNames = ['ContigHit','Start', 'End','Sequence','FIMOIscore']
    dfFullSeq = pd.DataFrame(ListSeq, columns = ColumnNames)
    dfFullSeq = dfFullSeq.astype({'FIMOIscore':'float'})
    return(dfFullSeq)

############################################
# TO CALCULATE GC CONTENT IN MOVING WINDOW #
############################################

def GCcalc(Sequence):
    Length = len(Sequence)
    GCcount = 0
    for letter in Sequence:
        if letter == 'G' or letter == 'C' or letter == 'c' or letter == 'g':
            GCcount = GCcount +1
    GCcont = GCcount/Length
    return(GCcont)

WindowSize = 100 #50-100 seems good

def GCwindow(Sequence,WindowSize):

    MovingGC = []
    Length = len(Sequence)

    for i in range(0,Length-WindowSize):
        Start = i
        End = i + WindowSize
        PlotCoord = End-(WindowSize/2)
        GCcont = GCcalc(Sequence[Start:End])
        newline = [Start, End, PlotCoord, GCcont]
        MovingGC.append(newline)

    GCdf = pd.DataFrame(MovingGC)
    GCdf.columns = ['Start','End','PlotCoord','GCcontent']
    return(GCdf)


################################################
# TO CALCULATE GC CONTENT OF CANDIDATE REGIONS #
################################################

dfFullSeq = FIMO2(Genome, Motif, Motifup, Motifdn)
DataFrames = pd.DataFrame()
os.mkdir("plots")

for index, row in dfFullSeq.iterrows():
    gcComp = GCwindow(row['Sequence'], WindowSize)
    dataFrame = gcComp
    
    fig = dataFrame.plot(x = "PlotCoord", y = "GCcontent", title = row['ContigHit'], ylim = (0,0.6))
    fig_title = row["ContigHit"]
    fig_title = fig_title[:fig_title.rfind("_")]
    fig.figure.savefig(f"plots/{fig_title}.png")
   
    dataFrame['HitNumber'] = index
    DataFrames = pd.concat([DataFrames, dataFrame])

