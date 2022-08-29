from nupack import *


def getPairProb_all(mySequence, my_model):
    pairsMatrix = pairs(strands=mySequence, model=my_model)
    pairsArray = pairsMatrix.to_array()
    return pairsArray

def getPairProb_individual(i, j, pairsArray):
    individualProb = pairsArray[i][j]
    return individualProb


