from nupack import *


def initializeRnaDict(numNucs):
    #make the dict and intialize each nuc pair to -1 or not paires. 
    #ensure that the smae nucs for i and j are retianed for investigations
    rnaDict={}
    for i_index in range(numNucs):
        rnaDict[i_index]={}
        for y_index in range(numNucs):
            rnaDict[i_index][y_index]=-1

    return rnaDict


def getPairProbsDict(mySequence, my_model):
    #get pairs array
    pairsDict =initializeRnaDict(len(mySequence))
    pairsMatrix = pairs(strands=mySequence, model=my_model)
    pairsArray = pairsMatrix.to_array()

    #now populate the pairs dict
    for i in range(len(pairsArray)):
        for j in range(len(pairsArray[i])):
            pairsDict[i][j]=pairsArray[i][j]

    return pairsDict

def getPairProb_individual(i, j, pairsDict):
    individualProb = pairsDict[i][j]
    return individualProb


