from varnaapi import *
from nupack import *


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