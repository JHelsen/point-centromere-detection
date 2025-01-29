#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 15:25:17 2024

@author: janahelsen

This script automates the retrieval of point centromere sequences from Saccharomycetaceae genome assemblies, using:
    - 2x FIMO from the MEME suite (Grant et al. 2011, Bioinformatics)
    - CDEII length and AT% based sorting and filtering
    - Synteny checks

"""

###########
# MODULES #
###########

import os
import re
import pandas as pd
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from orffinder import orffinder
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time
from tqdm import tqdm
from datetime import datetime

import warnings
warnings.filterwarnings('ignore', module='Bio')

####################
# BASIC PARAMETERS #
####################
PCAn_version = 0.1
print(f"\n\nWelcome to PCAn v{PCAn_version}, a tool to detect point centromeres in any genome assembly.")

#Input path to genome assembly file here
Genome = input("\nInput path to genome assembly file here. Please include the file name with extension in the path\n")#"/Users/janahelsen/Documents/Postdoc/Bioinformatics/Centromeres/Genomes/Zygotorulaspora_mrakii.fna"

print("\nNOTE: The tool uses the taxonomic classification as available on NCBI taxonomy prior to January 2025. Please ensure that you select the correct genus as per the current NCBI taxonomy.")
available_genus_list = ["Saccharomyces","Nakaseomyces","Oligophagozyma","Cylindricascospora","Arxiozyma","Monosporozyma","Grigorovia","Kazachstania","Huiozyma","Sungouiozyma","Maudiozyma","Jamesozyma","Naumovozyma","Tetrapisispora","Vanderwaltozyma","Yueomyces","Henningerozyma","Torulaspora","Zygotorulaspora","Zygosaccharomyces","Hagleromyces","Lachancea","Eremothecium","Ashbya","Kluyveromyces", "NA"]
for i in range(len(available_genus_list)):
    print(f"{i} - {available_genus_list[i]}")
del i
Species = available_genus_list[int(input("Please select the genus from the above:\n"))] #"Zygotorulaspora taianensis"
if Species == "Kazachstania":
    available_kazachstania_species = ['taianensis','bromeliacearum','martiniae','kunashirensis','psychrophila','viticola','africana']
    print("\n")
    for i in range(len(available_kazachstania_species)):        
        print(f"{i} - {available_kazachstania_species[i]}")
    del i
    print("\n")
    Species += " " + available_kazachstania_species[int(input("Please select the Kazachstania species from the above:\n"))]

if Species == "Yueomyces":
    print("WARNING: please note that for Yueomyces spp., a recognizable CDEI cannot be detected by PCAn. PCAn can give search results for the best-scoring CDEIII hits which can be found at CDEIII_search/SequenceHits\n(This folder is deleted for other species)\n")

if Species == "Naumovozyma":
    print(f"WARNING: Naumovozyma spp. have unconventional point centromeres. These spp. needed a different approach. In PCAn v{PCAn_version}, Nauvomozyma spp. centromeres cannot be automatically detected and require some further optimization.\n")

#Available options for genus names are Saccharomyces, Nakaseomyces, Oligophagozyma, Cylindricascospora, Arxiozyma, Monosporozyma, Grigorovia, Kazachstania, Huiozyma, Sungouiozyma, Maudiozyma, Jamesozyma, Naumovozyma, Tetrapisispora, Vanderwaltozyma, Yueomyces, Henningerozyma, Torulaspora, Zygotorulaspora, Zygosaccharomyces, Hagleromyces, Lachancea, Eremothecium, Ashbya, Kluyveromyces.
#These are the official genus names available at the time of publishing. In case of a new genus names, please use the most closely related genus available in the list above.
#If you have no idea about the phylogenetic placement, use "NA"

#!!! Should give warning if unknown

# Do you want to BLAST to check synteny? This will increase the computing time.
BLASTING = input("\nDo you want to BLAST to check synteny? [Yes/No]\n\nThis will increase the computing time\n")
print("\n")
#Input path to PCAn folder here
PCAnFolder = os.path.dirname(os.path.realpath(__file__))

#######################
# ADVANCED PARAMETERS #
#######################

# Dictionary with CDEI and CDEIII motif file names, and thresholds for FIMO searches specific for a genus or species.
#!!! Include option for custom motif
MotifsAndThresholds = {'Saccharomyces':["CDEIII/CDEIII_Saccharomyces_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Nakaseomyces':["CDEIII/CDEIII_Nakaseomyces_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Oligophagozyma':["CDEIII/CDEIII_ScSuCgVp_MEME.txt","CDEI/CDEI_ZT_MEME.txt",5,2],
                       'Cylindricascospora':["CDEIII/CDEIII_ScSuCgVp_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Arxiozyma':["CDEIII/CDEIII_Ktel_MEME.txt","CDEI/CDEI_Arxiozyma_MEME.txt",6,3],
                       'Monosporozyma':["CDEIII/CDEIII_Kaq-Ksol_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Grigorovia':["CDEIII/CDEIII_Grig_MEME.txt","CDEI/CDEI_Grig_MEME.txt",5,2],
                       'Kazachstania taianensis':["CDEIII/CDEIII_Ktai_MEME.txt","CDEI/CDEI_ZT_MEME.txt",7,2],
                       'Huiozyma':["CDEIII/CDEIII_Knagsin_MEME.txt","CDEI/CDEI_Knagsin_MEME.txt",5,2],
                       'Kazachstania bromeliacearum':["CDEIII/CDEIII_Ka_Full_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Sungouiozyma':["CDEIII/CDEIII_Ka_Full_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Kazachstania martiniae':["CDEIII/CDEIII_Ka_Full_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Maudiozyma':["CDEIII/CDEIII_Ka_Full_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Kazachstania kunashirensis':["CDEIII/CDEIII_Ka_Full_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Kazachstania psychrophila':["CDEIII/CDEIII_Ka_Full_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Jamesozyma':["CDEIII/CDEIII_Ka_Full_MEME.txt","CDEI/CDEI_Jamesozyma_MEME.txt",6,2],
                       'Kazachstania viticola':["CDEIII/CDEIII_Ka_Full_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Kazachstania africana':["CDEIII/CDEIII_Kafricana_MEME.txt","CDEI/CDEI_Kafricana_MEME.txt",6,2],
                       'Tetrapisispora':["CDEIII/CDEIII_Tetrapisispora_MEME.txt","CDEI/CDEI_Tetrapisispora_MEME.txt",5,2],
                       'Vanderwaltozyma':["CDEIII/CDEIII_Vanderwaltozyma_MEME.txt","CDEI/CDEI_Vanderwaltozyma_MEME.txt",5,2],
                       'Yueomyces':["CDEIII/CDEIII_Saccharomyces_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Henningerozyma':["CDEIII/CDEIII_Henningerozyma_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Torulaspora':["CDEIII/CDEIII_Toru_MEME.txt","CDEI/CDEI_Toru_MEME.txt",5,2],
                       'Zygotorulaspora':["CDEIII/CDEIII_ZygoToru_MEME.txt","CDEI/CDEI_ZT_MEME.txt",5,2],
                       'Zygosaccharomyces':["CDEIII/CDEIII_ZygoSacch_MEME.txt","CDEI/CDEI_ZygoSacch_MEME.txt",6,3],
                       'Hagleromyces':["CDEIII/CDEIII_KLE_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Lachancea':["CDEIII/CDEIII_Lachancea_MEME.txt","CDEI/CDEI_Lachancea_MEME.txt",5,2],
                       'Eremothecium':["CDEIII/CDEIII_EreLach_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Ashbya':["CDEIII/CDEIII_EreLach_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2],
                       'Kluyveromyces':["CDEIII/CDEIII_Kluyveromyces_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,3],
                       'NA':["CDEIII/CDEIII_ScSuCgVp_MEME.txt","CDEI/CDEI_ZT_MEME.txt",6,2]}

for i in MotifsAndThresholds:
    MotifsAndThresholds[i][0] = os.path.join(PCAnFolder, MotifsAndThresholds[i][0])
    MotifsAndThresholds[i][1] = os.path.join(PCAnFolder, MotifsAndThresholds[i][1])
del i

# How many base pairs upstream of CDEIII to look for CDEI
CDEIIIup = 250 #default = 250

# How many top-ranked hits to retain in the first filtering step
Nrows = 50 #default = 50

# How many base pairs up- and downstream of potential centromere to extract for for synteny BLASTing
Range = 10000 #default = 15000

# Minimum number of amino acids for ORF prediction
MinORFLength = 175 #default = 175

################
# INITIALIZING #
################ 

# Starting the clock
startTime = time.localtime()

# Initializing report
Report = []

# Setting the working directory
# os.chdir(PCAnFolder)

####################################
# DEFINING THE TWICE FIMO FUNCTION #
####################################

def FIMO2(Genome, CDEI, CDEIII, CDEIIIup, TH1, TH2):

    command1 = "fimo --thresh 1.0E-"+ str(TH1) +" --oc FIMOresults "+ CDEIII + " "+ Genome
    os.system(command1)
    
    fimosearchResult = 'FIMOresults/fimo.txt'
    Queries = pd.read_csv(fimosearchResult, skiprows=0,delimiter="\t")
    
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
            Newstart = start - CDEIIIup
            Newend = end
            POSITION1 = Newend
        if strand == '-':
            Newstart = start -1
            Newend = end + CDEIIIup
            POSITION1 = Newstart
        for seq in SeqIO.parse(Genome,"fasta"):
            if seq.id == ID:
                if strand == '+':
                    sequence = str(seq.seq[Newstart:Newend])
                if strand == '-':
                    seqRC = seq.seq[Newstart:Newend]
                    sequence = str(seqRC.reverse_complement())          
                Sequences.append(sequence)
        Lines.append(">"+ ID + "_" + str(index+1) + "\n" + sequence)
        ListSeq.append([ID + "_" + str(index+1), sequence, QualityScoreI,POSITION1,strand])
    
    os.mkdir("CDEIII_search")
    
    with open('CDEIII_search/SequenceHits.fsa','a') as outfile:
        outfile.writelines(i + '\n' for i in Lines)
        outfile.close
    
    command2 = "fimo --norc --thresh 1.0E-"+ str(TH2) +" --oc FIMOresults " + CDEI + " "+os.getcwd()+"/CDEIII_search/SequenceHits.fsa"
    os.system(command2)
    
    if Species != "Yueomyces":
        os.remove("CDEIII_search/SequenceHits.fsa")
        os.rmdir("CDEIII_search")
    
    with open(fimosearchResult, 'r') as M:
        reader = csv.reader(M, delimiter = '\t')
        next(reader)
        FMhits = list(reader)
    
    FullSeq = []
    
    for FMhit in FMhits:
        Start = int(FMhit[2])
        End = int(FMhit[3])
        for Hit in ListSeq:
            if FMhit[1] == Hit[0]:
                Contig = re.search('.+(?=_.+)',FMhit[1]).group(0)
                ContigHit = Hit[0]
                ShortSeq = Hit[1][Start-1:]
                ATSeq = Hit[1][End+1:len(Hit[1])-25]
                ATlen = len(ATSeq)
                ATper = round((ATSeq.count("A")+ATSeq.count("a")+ATSeq.count("T")+ATSeq.count("t"))/(ATlen+1)*100, 2)
                QualityScoreI = Hit[2]
                QualityScoreII = float(FMhit[5])
                OverallQuality = round(QualityScoreI*QualityScoreII*1/(1-ATper/100), 0)
                if Hit[4] == "+":
                    NewStart = Hit[3] - len(ShortSeq) + 1
                    NewEnd = Hit[3]
                else:
                    NewStart = Hit[3] + 1
                    NewEnd = Hit[3] + len(ShortSeq)
                newline = [Contig, ContigHit, NewStart, NewEnd, ShortSeq , str(ATlen) , str(ATper), str(QualityScoreI), str(QualityScoreII), str(OverallQuality)]
                FullSeq.append(newline)
    
    ColumnNames = ['Contig', 'ContigHit', 'Start', 'End','Sequence','CDEIIlen','CDEIIAT','FIMOIscore','FIMOIIscore','OverallScore']
    dfFullSeq = pd.DataFrame(FullSeq, columns = ColumnNames)
    dfFullSeq = dfFullSeq.astype({'Start':'int','End':'int','CDEIIlen':'int','CDEIIAT':'float','FIMOIscore':'float','FIMOIIscore':'float','OverallScore':'float'})
    SortedFull = dfFullSeq.sort_values(by='OverallScore',ascending = False, ignore_index = True)
    return(SortedFull)

######################
# RUNNING TWICE FIMO #
######################

# Chosing parameters based on species
ChosenParameters = MotifsAndThresholds[Species]

dfFullSeq = FIMO2(Genome, ChosenParameters[1], ChosenParameters[0], CDEIIIup, ChosenParameters[2], ChosenParameters[3])   

AFterFIMOTime = time.localtime()
AnalysisTime = (time.mktime(AFterFIMOTime) - time.mktime(startTime))/60
Atime = "The total FIMO analysis time was " + str(AnalysisTime) + " minutes."
Report.append(Atime)


############################################
# FILTERING BASED ON LENGTH AND SUB-CONTIG #
############################################

def Filtering(CENdf):

    MedianCDEIIlen = dfFullSeq.head(5)['CDEIIlen'].median()
    
    # Selecting top
    dfHead = dfFullSeq.head(Nrows)
    
    # Selecting all hits around median
    dfHeadCDEIIlen = dfHead.loc[(dfHead['CDEIIlen'] > MedianCDEIIlen - 30) & (dfHead['CDEIIlen'] < MedianCDEIIlen + 30)]
    
    # Removing low AT%
    
    dfHeadAT = dfHeadCDEIIlen.loc[dfHeadCDEIIlen['CDEIIAT'] > 70]

    # Removing lines with negative motif scores
  
    dfMotifPos = dfHeadAT.loc[(dfHeadAT['FIMOIscore'] > 0) & (dfHeadAT['FIMOIIscore'] > 0)]

    # Removing doublets: same CDEIII hit, but different CDEI hit
    
    dfHeadCDEIIlenDoub = dfMotifPos.drop_duplicates(subset = 'ContigHit', keep = 'first')
    
    # Removing doublets: same exact CEN
    
    dfHeadCENDoub = dfHeadCDEIIlenDoub.drop_duplicates(subset = 'Sequence', keep = 'first')
    
    # Removing doublets: same contig, multiple hits
    
    dfHeadMult = dfHeadCENDoub.drop_duplicates(subset = 'Contig', keep = 'first')
        
    # Removing lines with outliers in both AT% and len
    
    Medlen = dfHeadMult['CDEIIlen'].median()
    MedAT = dfHeadMult['CDEIIAT'].quantile(0.5)
    
 #!!! Length used to be +/-5
    dfHeadAllFilt = dfHeadMult.loc[((dfHeadMult['CDEIIlen'] > Medlen - 10) & (dfHeadMult['CDEIIlen'] < Medlen + 10)) | (dfHeadMult['CDEIIAT'] > MedAT - 7)]
    
    return(dfHeadAllFilt)

FilteredDf = Filtering(dfFullSeq)

CENNo = str(len(FilteredDf)) + " potential centromeres were identified."
Report.append(CENNo)

OutputFileName = os.getcwd()+"/CENsequences.txt"

FilteredDf.to_csv(OutputFileName, sep = '\t', index = False)

with open(OutputFileName, 'a') as f:
    f.write("\n****************\n\n")
    f.write(f"PCAn v{PCAn_version}\n")
    f.write(datetime.now().strftime('%d-%m-%Y %H:%M:%S'))
    f.write(f"\n\n{CENNo}\n")
    f.write(f"\nPath to genome used: {Genome}\n")
    f.write(f"Genus/species selected: {Species}\n")
    f.write(f"\nCDEI motif used: {MotifsAndThresholds[Species][1]}\n")
    f.write(f"FIMO significance threshold for the CDEI motif search: 1E-{MotifsAndThresholds[Species][3]}\n")
    f.write(f"\nCDEIII motif used: {MotifsAndThresholds[Species][0]}\n")
    f.write(f"FIMO significance threshold for the CDEIII motif search: 1E-{MotifsAndThresholds[Species][2]}\n")
    f.write("\nPlease cite:\nJ. Helsen, K. Ramachandran, G. Sherlock, G. Dey, Centromeres evolve progressively through selection at the kinetochore interface. bioRxiv [Preprint] (2025). https://doi.org/10.1101/2025.01.16.633479.")
    f.close()

###################################
# BLAST SEARCHES ON FINAL RESULTS #
###################################

#!!! BLASTin and BLASTout folders are not emptied, made or removed

# Function that returns the sequence x bp up and down (see range), given a CEN sequence, contig, and genome
# The function also returns the location of the CEN sequence within the new extracted sequence

if BLASTING == 'Yes':
    print("\n\nRunning BLASTp searches for synteny analysis\n")
# Extracting neighbouring sequences    

    def extr(seq,cont,genome,updn):
        for C in SeqIO.parse(genome,"fasta"):
            C_id = C.id
            if C_id == cont:
                Sequence = C.seq.upper()
                CompSequence = C.reverse_complement().seq.upper()
        if Sequence.find(seq.upper()) != -1:
            Coor = Sequence.find(seq.upper())
            Correctseq = Sequence
        else:
            Coor = CompSequence.find(seq.upper())
            Correctseq = CompSequence
        StartUp = Coor - 1 - updn
        if StartUp < 0:
            StartUp = 0
            StartCEN = StartUp + updn
        else:
            StartCEN = updn
        EndDn = Coor + len(seq) + 1 + updn
        if EndDn > len(Correctseq):
            EndDn = -1
        EndCEN = StartCEN + len(seq)
        WHOLE = Correctseq[StartUp:EndDn]
        record = SeqRecord(Seq(str(WHOLE)))
        return(record,[StartCEN,EndCEN])
    
    # Looping over all CEN sequences in the list
    
    Synteny = []
    
    for ind in FilteredDf.index:
        contig = FilteredDf['Contig'][ind]
        x = extr(FilteredDf['Sequence'][ind], contig, Genome, Range)
        newline = []
        newline.append(contig)  # Contig
        newline.append(x[0])    # Whole seq
        newline.append(x[1])    # CEN position in seq
        Synteny.append(newline)  
        
    # Finding ORFs and BLASTing
    
    # Making the database (needs to be done once)
    #cline = NcbimakeblastdbCommandline(dbtype="prot", input_file="BLASTdb/AllScORFs.fasta", cmd= '/usr/local/ncbi/blast/bin/makeblastdb')
    #cline()
    
    # Short function that recognizes empty files in case of no blast results
    def is_non_zero_file(fpath):  
        return os.path.isfile(fpath) and os.path.getsize(fpath) > 0 
    
    os.mkdir("BLAST")
    os.chdir("BLAST")
    os.mkdir("BLASTin")
    os.mkdir("BLASTout")
    os.chdir(os.path.dirname(os.getcwd()))

    for CEN in tqdm(Synteny, unit=" centromeres", position=0, leave=True):
        ORFs = orffinder.getORFProteins(CEN[1], minimum_length=MinORFLength*3, remove_nested= True, return_loci=True)
        FoundORFs = []
        for i in range(0,len(ORFs)):
            Name = CEN[0] + '.' + str(i)
            Sequence = ORFs[i]['protein']
            record = SeqRecord(Sequence,id=Name)
            PathIn = 'BLAST/BLASTin/' + Name + '.fsa'
            PathOut = 'BLAST/BLASTout/' + Name + '.txt'
            SeqIO.write(record, PathIn,'fasta')
            # cmd_blastp = NcbiblastpCommandline(query = PathIn, out = PathOut, outfmt = 6, db = os.path.join(PCAnFolder, "cerevisiae_proteome/Sc_proteome.fasta"), cmd = 'blastp')
            # cmd_blastp()
            os.system(f'blastp -db {PCAnFolder + "/cerevisiae_proteome/Sc_Proteome.fasta"} -query {PathIn} -out {PathOut} -outfmt 6 >/dev/null 2>&1')
            if is_non_zero_file(PathOut):
                results = pd.read_csv(PathOut, sep="\t", header=None)
            else:
                results = [['-'],['NoHit']]
            entry = [Name,ORFs[i]['start'],ORFs[i]['end'],results[1][0]]
            FoundORFs.append(entry)
            # os.remove(PathIn)
            # os.remove(PathOut)
        CEN.append(FoundORFs)
    print("BLASTp searches complete\n")    
    # Rearranging and adding ancestral protein IDs
    
    OrderedSyn = []
    
    for CEN in Synteny:
        newline = []
        newline.append(CEN[0])
        centers = []
        CENPOS = [CEN[0],CEN[2][0],CEN[2][1],'CEN']
        withCEN = [i for i in CEN[3]]
        withCEN.append(CENPOS)
        for ORF in withCEN:
            center = (ORF[1]+ORF[2])/2
            centers.append(center)
        sorted_Ind = [i[0] for i in sorted(enumerate(centers), key=lambda x:x[1])]
        sorted_ORFs = []
        for i in sorted_Ind:
            sorted_ORFs.append(withCEN[i][3])
        newline.append(sorted_ORFs)
        OrderedSyn.append(newline)
     
        
    with open(os.path.join(PCAnFolder, 'ANCPIL.txt'), 'r') as M:
        reader = csv.reader(M, delimiter = '\t')
        ANCPIL = list(reader)
    
        
    for CEN in OrderedSyn:
        ancORFs = []
        for hit in CEN[1]:
            if hit != 'CEN' and hit!= 'NoHit':
                for line in ANCPIL:
                    if hit == line[-1] or hit == line[-2]:
                        ancORFs.append(line[1])
            elif hit == 'CEN':
                ancORFs.append('CEN')
            else:
                ancORFs.append('NoHit')
        CEN.append(ancORFs)
    
    # Plotting output
    
    fig, ax = plt.subplots()
    
    for j in range(0,len(OrderedSyn)):
        ax.text(-0.3, j*0.15, OrderedSyn[j][0], color = 'black', bbox=dict(facecolor='none', edgecolor='black', boxstyle = 'square'))
        for i in range(0,len(OrderedSyn[j][2])):
            if OrderedSyn[j][2][i] == 'CEN':
                C = 'black'
                FC = 'white'
                BS = 'circle'
            if OrderedSyn[j][2][i] == 'NoHit':
                C = 'black'
                FC = 'white'
                BS = 'square'
            elif re.search('._L.+', OrderedSyn[j][2][i]):
                C = 'blue'
                FC = 'dodgerblue'
                BS = 'square'
            elif re.search('._R.+', OrderedSyn[j][2][i]):
                C = 'darkmagenta'
                FC = 'magenta'
                BS = 'square'
            ax.text(i*0.15, j*0.15, OrderedSyn[j][2][i], color = 'black', bbox=dict(facecolor=FC, edgecolor=C, boxstyle = BS))
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.xaxis.set_major_locator(ticker.NullLocator())
            ax.yaxis.set_major_locator(ticker.NullLocator())
    
    plt.savefig('SyntenyPlot.png', bbox_inches = 'tight', dpi = 160)

endTime = time.localtime()
AnalysisTime = (time.mktime(endTime) - time.mktime(startTime))/60
Atime = "The total analysis time was " + str(AnalysisTime) + " minutes."
Report.append(Atime)


ReportFileName = os.path.join(os.getcwd(), "Report.txt")


