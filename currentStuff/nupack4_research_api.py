from nupack import *

my_model = Model

rna_model='rna95-nupack3'    
# Define physical model 



def SetModel(rna_model, temp_C):
    my_model = Model(material=rna_model, celsius=int(temp_C))
    return my_model

def initializeRnaDict(numNucs):
    #make the dict and intialize each nuc pair to -1 or not paires. 
    #ensure that the smae nucs for i and j are retianed for investigations
    rnaDict={}
    for i_index in range(numNucs):
        rnaDict[i_index]={}
        for y_index in range(numNucs):
            rnaDict[i_index][y_index]=-1

    return rnaDict


def getPairProbs(mySequence, temp,  filterCutoff):
    #get pairs array
    theModel = SetModel(rna_model, temp)
    pairsDict =initializeRnaDict(len(mySequence))
    pairsMatrix = pairs(strands=mySequence, model=theModel)
    pairsArray = pairsMatrix.to_array()
    newPairsList=list()
    snupp_PaisList = []
    #now populate the pairs dict
    for i in range(len(pairsArray)):
        newPairsList.append(list())
        for j in range(len(pairsArray[i])):
            value = pairsArray[i][j]
            if filterCutoff=="l2":
                if  value < .01 and value >= .001:
                    if i is not j:
                        pairsDict[i][j]=value
                        newPairsList[i].append(value)
                        snupPairValue = "{0}:{1}={2:.10f}".format(i+1, j+1, value)
                        snupp_PaisList.append(snupPairValue)
            else:
                if  value >= float(filterCutoff):
                    if i is not j:
                        pairsDict[i][j]=value
                        newPairsList[i].append(value)
                        snupPairValue = "{0}:{1}={2:.10f}".format(i+1, j+1, value)
                        snupp_PaisList.append(snupPairValue)

    return pairsDict, newPairsList, snupp_PaisList

def getPairProb_individual(i, j, pairsDict):
    individualProb = pairsDict[i][j]
    return individualProb


#structure stuff



