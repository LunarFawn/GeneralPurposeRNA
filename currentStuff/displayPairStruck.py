#from varnaapi import *
import struct
from nupack import *
import pandas as pd
import sys
import openpyxl


print("Enter single strand RNA sequence")
sequence = "GGGAACGACUCGAGUAGAGUCGAAAAGUUGAAACGACUGAACAUGGUAACAUGAGUGGUUUGAACACAUACGAACAGGGUUCUUUCGAGGAUCCAAAAGAAACAACAACAACAAC" #input()

print("Enter the temperature to simulate at")
temperature = 37 #input()

rna_model='rna95-nupack3'    
# Define physical model 

def energy_model(rna_model, temperature):
    model = Model(material=rna_model, celsius=int(temperature))
    return model

my_model = energy_model(rna_model, temperature)

def processPairs(mySequence, primaryCutoff, secondaryCuttoff):
    
    pairsMatrix = pairs(strands=mySequence, model=my_model)
    pairsArray = pairsMatrix.to_array();
    #make a two dimensional list that has all teh pairs for everything
    nucs = len(mySequence)
    nucPairsList=[[-1]*nucs]*nucs
        
    #now go through the probability matrix and placie it in the right place
    primaryPairsSortedList=[]
    primaryPairsList=[]
    print("")
    print("Primary Pairs List")
    for i in range(len(pairsArray)):
        for j in range(len(pairsArray[i])):
            if pairsArray[i][j]>primaryCutoff:
                message = "{0:0>2d}:{1:0>2d}={2:.6f}".format(i+1, j+1, pairsArray[i][j])
                primaryPairsList.append(message)
                sortmessage = "{2:.6f}={0:0>2d}:{1:0>2d}".format(i+1, j+1, pairsArray[i][j])
                #bisect.insort(primaryPairsList, sortmessage)
                primaryPairsSortedList.append(sortmessage)
                print(message)
    print("")
    print("Sorted Primary Pairs")
    primaryPairsSortedList.sort(reverse=True)
    for element in primaryPairsSortedList:
        print(element)
    
    secondaryPairsSortedList=[]
    secondaryPairsList=[]                
    print("")
    print("Secondary Pairs List")
    for i in range(len(pairsArray)):
        for j in range(len(pairsArray[i])):
            if pairsArray[i][j]<primaryCutoff and pairsArray[i][j]>secondaryCuttoff:
                message = "{0:0>2d}:{1:0>2d}={2:.6f}".format(i+1, j+1, pairsArray[i][j])
                secondaryPairsList.append(message)
                sortmessage = "{2:.6f}={0:0>2d}:{1:0>2d}".format(i+1, j+1, pairsArray[i][j])
                #bisect.insort(secondaryPairsList, sortmessage)
                secondaryPairsSortedList.append(sortmessage)
                print(message)
    print("")
    print("Sorted Secondary Pairs")
    secondaryPairsSortedList.sort(reverse=True)
    for element in secondaryPairsSortedList:
        print(element)
    #now return all the lists for ingesting
    return pairsArray, primaryPairsList, primaryPairsSortedList, secondaryPairsList, secondaryPairsSortedList



#now check for 100% and also specified fold change
def grabPnasData():
    pnasPath = r'pnas.2112979119.sd01.xlsx'
    designRound = "Round 7 (R101)"
    
    sheet = pd.read_excel(pnasPath, sheet_name=designRound)

    index=1
    for lab_design_index in sheet.index:
        sequence = sheet['Sequence'][lab_design_index]
        DesignID = sheet['DesignID'][lab_design_index]
        DesignName = sheet['Design'][lab_design_index]
        Player = sheet['Player'][lab_design_index]
        Eterna_Score=sheet['Eterna_Score'][lab_design_index]
        FoldChange=sheet['FoldChange'][lab_design_index]
        Puzzle_Name=sheet['Puzzle_Name'][lab_design_index]
        #run single sequence
        
       



        test = new_dict.get_entry('sequence')
        result = test.get_value()
