#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 08:41:24 2024

@author: janahelsen
"""

##############################
# Simulating CEN transitions #
##############################

import numpy as np
rng = np.random.default_rng(19680801) #number is the seed
import pandas as pd

########################
# HAPLOIDS WITHOUT SEX #
########################

Ntransitions = 10
Nsimulations = 10000

StartState = ['A']*16


###    
# Defining function for 1 simulation - forcing state transition, but varying the chance of retention within the population
###

def simfitness(startState, noTransitions, threshold):
    MutationRecord = [startState]
    CENs = []
    for i in range(0,noTransitions):
        newState = []
    # Picking a random CEN     
        CEN = rng.integers(16,size=1)
    # Fitness effect of mutation
        MUT = rng.integers(1000, size = 1)    
        for j in range(0,len(startState)):
            if j != CEN[0]:
                newState.append(MutationRecord[i][j])
            else:
                oldCEN = MutationRecord[i][j]
                if oldCEN == 'A':
                    if MUT[0] <= threshold: #Here I'm implementing that A is better, so to change to B you need to be over the threshold
                        newCEN = 'A'
                    else:
                        newCEN = 'B'
                if oldCEN == 'B':
                    if MUT[0] > threshold: #Here I'm implementing that B is worse, so to keep B you need to be over the threshold
                        newCEN = 'B'
                    else:
                        newCEN = 'A'
                newState.append(newCEN)
        MutationRecord.append(newState)
        CENs.append(CEN[0])
    return(MutationRecord[-1])
        
# Doing many simulations

SimulationRecord = []

for i in range(0,Nsimulations):
    SIM = simfitness(StartState,Ntransitions,999)
    CENcount = SIM.count('A')
    SimulationRecord.append(CENcount)

CompiledRecord = []

for i in range (0,17):
    count = SimulationRecord.count(i)
    CompiledRecord.append(count)
    

###############################################
# DIPLOIDS WITH SEX - SET NUMBER OF MUTATIONS #
###############################################

NMut = 500 # This is a mutation event
Nsimulations = 2000
PopulationSize = 100 # Size will be kept constant.

StartState = [['A', 'A']]*16
StartPop = [StartState]*PopulationSize

MutRetentionThreshold = 0 #Threshold for mutation retention

PopThreshold = 100 #Number of individuals within the population (depends on population size).


def simFitnessSex(startPop, nTransitions, mutRetThres, popThres, recombination):
    MutationRecord = startPop
    # Looping over all transitions
    for i in range(0,nTransitions):       
        newPop = []
        # Choosing random individuals within the population by shuffling indices
        indices = list(range(len(startPop)))
        rng.shuffle(indices)
        Chosen = indices[:popThres]
        # Looping over all individuals
        for k in range(0, len(startPop)):
            OldInd = MutationRecord[k]
            # If not chosen for a transition
            if k not in Chosen:
                NewInd = OldInd
            # If chosen individual
            if k in Chosen:
                NewInd = []
                # Picking a random CEN     
                CEN = rng.integers(16, size = 1)
                # Picking a random chromatid
                CHD = rng.integers(2, size = 1)
                # Fitness effect of mutation
                MUT = rng.integers(100, size = 1)   
                # Looping over centromeres
                for j in range(0,len(OldInd)):
                    # If centromere wasn't chosen: append old variant
                    if j != CEN[0]:
                        NewInd.append(OldInd[j])
                    # If centromere was chosen
                    else:
                        newCEN_xy = []
                        oldCEN_xy = OldInd[j]
                        for chrtd in [0,1]:
                            # If chosen chromatid: mutate
                            if chrtd == CHD[0]:
                                if oldCEN_xy[chrtd] == 'A':
                                    if MUT[0] <= mutRetThres: #Here I'm implementing that A is better, so to change to B you need to be over the threshold
                                        newCEN_x = 'A'
                                    else:
                                        newCEN_x = 'B'
                                if oldCEN_xy[chrtd] == 'B':
                                    if MUT[0] > mutRetThres: #Here I'm implementing that B is worse, so to keep B you need to be over the threshold
                                        newCEN_x = 'B'
                                    else:
                                        newCEN_x = 'A'
                                newCEN_xy.append(newCEN_x)
                            # If not chosen chromatid: use old variant
                            else:
                                newCEN_xy.append(oldCEN_xy[chrtd])
                        # Append new variant
                        NewInd.append(newCEN_xy)
            if recombination == 'Yes':
                # Recombination
                AfterRec = []
                # For each CEN
                for cen in range(0, len(NewInd)):
                    newCEN_xy = []
                    oldCEN_xy = NewInd[cen]
                    # Picking a random chromatid to retain
                    CHDrec = rng.integers(2, size = 1)
                    for chromatid in [0,1]:
                        # If chosen chromatid: sex - choose a random variant from the population 
                        if chromatid == CHDrec[0]:
                            chosenExchangePartner = rng.integers(len(MutationRecord), size = 1)
                            chosenExchangeChrom = rng.integers(2, size = 1)
                            newCEN_x = MutationRecord[chosenExchangePartner[0]][cen][chosenExchangeChrom[0]]
                        else:
                            newCEN_x = oldCEN_xy[chromatid]
                        newCEN_xy.append(newCEN_x)
                    AfterRec.append(newCEN_xy)
                newPop.append(AfterRec)
            else:
                newPop.append(NewInd)
        MutationRecord = newPop         
    return(MutationRecord)
    
# Doing many simulations

SimulationRecordA = []
SimulationRecordB = []
IndividualAllA = [['A', 'A']]*16
IndividualAllB = [['B', 'B']]*16

for i in range(0,Nsimulations):
    SIM = simFitnessSex(StartPop,NMut,MutRetentionThreshold,PopThreshold, 'No')
    # Counting how many times all A or all B
    PopcountA = SIM.count(IndividualAllA)
    PopcountB = SIM.count(IndividualAllB)
    SimulationRecordA.append(PopcountA)
    SimulationRecordB.append(PopcountB)

AAcompiled = [[x, SimulationRecordA.count(x)] for x in set(SimulationRecordA)]
BBcompiled = [[x, SimulationRecordB.count(x)] for x in set(SimulationRecordB)]



AA = pd.DataFrame(AAcompiled, columns = ['NoIndPop', 'NoSim'])
BB = pd.DataFrame(BBcompiled, columns = ['NoIndPop', 'NoSim'])



